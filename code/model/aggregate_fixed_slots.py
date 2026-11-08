"""Downsample published interval speeds to complete 15-minute pace snapshots.

Each source interval is piecewise constant. Overlap-weighted arithmetic pace
means a harmonic mean of speed. Never upsample a coarser source interval.
"""
import argparse,csv,datetime,hashlib,json,shutil
from pathlib import Path
import numpy as np

DT=np.dtype([('e','<i4'),('w','<i4'),('kph','<f4')])

def overlap_plan(tags,source_seconds,target_seconds=900):
    if source_seconds>target_seconds:raise ValueError('Upsampling is forbidden: source intervals exceed target slots')
    tags=np.asarray(tags,dtype=np.int64)
    if len(np.unique(tags))!=len(tags):raise ValueError('Duplicate source slot tags')
    plans=[]
    for slot in range(86400//target_seconds):
        start=slot*target_seconds;stop=start+target_seconds
        weights=np.maximum(0,np.minimum(tags*source_seconds+source_seconds,stop)-np.maximum(tags*source_seconds,start))
        take=np.flatnonzero(weights)
        if int(weights.sum())==target_seconds:plans.append((slot,take,weights[take]))
    return plans

def aggregate_pace(speeds,seconds,target_seconds=900):
    values=np.asarray(speeds,dtype=np.float64);weights=np.asarray(seconds,dtype=np.float64)[:,None]
    valid=np.isfinite(values)&(values>0)
    inverse=np.divide(1.,values,out=np.zeros_like(values),where=valid)
    coverage=(valid*weights).sum(axis=0)
    inv_sum=(inverse*weights).sum(axis=0)
    out=np.full(values.shape[1],np.nan,dtype=np.float32)
    good=(coverage==target_seconds)&(inv_sum>0)
    out[good]=(target_seconds/inv_sum[good]).astype(np.float32)
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--out',required=True)
    p.add_argument('--seed',type=int,default=20261105);p.add_argument('--holdout-fraction',type=float,default=.05);a=p.parse_args()
    src=Path(a.input);out=Path(a.out)
    if out.exists():raise FileExistsError(out)
    n,m,days,slots,source_seconds,lo,hi,nnz,classes=map(int,(src/'meta.txt').read_text().split())
    count=hi-lo;oldN=days*count
    tags=np.fromfile(src/'slots.bin',dtype='<i4') if (src/'slots.bin').exists() else np.arange(lo,hi,dtype=np.int32)
    assert len(tags)==count and len(set(tags))==count
    plan=overlap_plan(tags,source_seconds)
    if not plan:raise ValueError('No fully covered 15-minute slots')
    obs=np.memmap(src/'obs.bin',dtype=DT,mode='r',shape=(nnz,))
    assert np.all((obs['e']>=0)&(obs['e']<m)&(obs['w']>=0)&(obs['w']<oldN))
    assert np.all(np.diff(obs['w'])>=0),'Source observations must be grouped by snapshot'
    speed=np.full((oldN,m),np.nan,dtype=np.float32)
    offsets=np.r_[0,np.cumsum(np.bincount(obs['w'],minlength=oldN))]
    for w in range(oldN):
        block=obs[offsets[w]:offsets[w+1]]
        if len(np.unique(block['e']))!=len(block):raise ValueError('Duplicate edge/source-snapshot observation')
        speed[w,block['e']]=block['kph']
    newN=days*len(plan);aggregate=np.empty((newN,m),dtype=np.float32);mapping=[]
    for day in range(days):
        for j,(slot,positions,weights) in enumerate(plan):
            target=day*len(plan)+j;source_worlds=day*count+positions
            aggregate[target]=aggregate_pace(speed[source_worlds],weights)
            mapping.append([target,day,slot,';'.join(map(str,source_worlds)),';'.join(map(str,weights)),
                            int(np.count_nonzero(np.isfinite(aggregate[target])))])
    valid=np.isfinite(aggregate)&(aggregate>0)
    ww,ee=np.nonzero(valid)
    # Hold out entire edge-day blocks, so overlapping source intervals cannot
    # occur in both training and held-out quarter-hour entries of the same edge.
    rng=np.random.default_rng(a.seed);blocks=rng.random((days,m))<a.holdout_fraction
    held=blocks[ww//len(plan),ee].astype(np.uint8)
    train=valid.copy();train[ww[held==1],ee[held==1]]=False
    ff=np.full(m,np.nan,dtype=np.float32)
    for e in range(m):
        if train[:,e].any():ff[e]=np.quantile(aggregate[train[:,e],e],.95)
    fallback=int(np.count_nonzero(~np.isfinite(ff)))
    if not np.isfinite(ff).any():raise ValueError('No training speeds for free-flow estimation')
    ff[~np.isfinite(ff)]=np.nanmedian(ff)
    out.mkdir(parents=True)
    records=np.empty(len(ee),dtype=DT);records['e']=ee;records['w']=ww;records['kph']=aggregate[ww,ee]
    records.tofile(out/'obs.bin');held.tofile(out/'holdout.bin');blocks.astype(np.uint8).tofile(out/'holdout_edge_day.bin')
    arcs=(src/'arcs.bin').read_bytes();assert len(arcs)==20*m
    (out/'arcs.bin').write_bytes(arcs[:12*m]+ff.astype('<f4').tobytes()+arcs[16*m:])
    for name in ['coords.bin','days.bin']:shutil.copyfile(src/name,out/name)
    np.array([x[0] for x in plan],dtype='<i4').tofile(out/'slots.bin')
    (out/'meta.txt').write_text(f'{n} {m} {days} 96 900 0 {len(plan)} {len(ee)} {classes}\n')
    if (src/'edge_mapping.csv').exists():
        with (src/'edge_mapping.csv').open() as f:
            reader=csv.DictReader(f);columns=reader.fieldnames;rows=list(reader)
        for e,row in enumerate(rows):
            if 'estimated_freeflow_kph' in row:row['estimated_freeflow_kph']=str(float(ff[e]))
        with (out/'edge_mapping.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
    with (out/'aggregation_map.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['target_world','day_index','slot_15min','source_worlds','overlap_seconds','valid_edges']);w.writerows(mapping)
    provenance=[]
    for name in ['meta.txt','obs.bin','arcs.bin','slots.bin','days.bin','preparation.json']:
        file=src/name
        if file.exists():
            with file.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
            provenance.append(dict(path=str(file),bytes=file.stat().st_size,sha256=sha))
    info=json.loads((src/'preparation.json').read_text())
    info.update(source_input_directory=str(src),source_slot_seconds=source_seconds,source_grid_worlds=oldN,
        N=newN,slot_seconds=900,slots_per_day=96,observations=len(ee),density=float(valid.mean()),
        heldout=int(held.sum()),heldout_fraction=float(held.mean()),empty_worlds=int(np.count_nonzero(valid.sum(axis=1)==0)),
        aggregation='sum(overlap_seconds * source_pace) / 900; each published source interval treated as constant',
        missing_policy='require all 900 seconds for an edge-slot; no interpolation across missing intervals',
        holdout_unit='edge-day block; all quarter-hour entries of that edge-day share one assignment',
        holdout_seed=a.seed,requested_holdout_fraction=a.holdout_fraction,
        freeflow='95th percentile of aggregated training speeds; held-out edge-day blocks excluded',
        freeflow_fallback_edges=fallback,one_snapshot_per_day_slot=True,upsampling=False,
        source_files_hashed=provenance,
        temporal_scope='Derived 15-minute averages of published interval estimates, not new GPS observations',
        min_positive_kph=float(aggregate[valid].min()),max_positive_kph=float(aggregate[valid].max()))
    (out/'preparation.json').write_text(json.dumps(info,indent=2))
    print(json.dumps({k:info[k] for k in ['city','n','m','days','N','slot_seconds','observations','heldout_fraction','empty_worlds','freeflow_fallback_edges']}))

if __name__=='__main__':main()
