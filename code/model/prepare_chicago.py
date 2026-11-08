"""Prepare official irregular Chicago speed reports under the frozen 15-minute rules."""
import argparse,csv,datetime,hashlib,json,math,re
from collections import defaultdict
from pathlib import Path
import numpy as np

def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def aggregate_reports(times,speeds):
    times=np.asarray(times,dtype=float);speeds=np.asarray(speeds,dtype=float)
    if len(times)!=len(speeds) or np.any(np.diff(times)<=0):raise ValueError('reports must have unique increasing times')
    if len(times) and (times[0]<0 or times[-1]>=86400):raise ValueError('report outside its calendar date')
    coverage=np.zeros(96);pace=np.zeros(96)
    for i,(start,speed) in enumerate(zip(times,speeds)):
        stop=min(start+900,times[i+1] if i+1<len(times) else 86400,86400)
        if not math.isfinite(speed) or speed<=0:continue
        for slot in range(max(1,int(start//900)),min(95,int((stop-1e-8)//900))+1):
            seconds=max(0,min(stop,(slot+1)*900)-max(start,slot*900))
            coverage[slot]+=seconds;pace[slot]+=seconds/speed
    valid=(np.abs(coverage-900)<1e-7)&(pace>0);out=np.full(96,np.nan,dtype=np.float32)
    out[valid]=(900/pace[valid]*1.609344).astype(np.float32)
    return out

def largest_scc(n,src,dst):
    adj=[[] for _ in range(n)];rev=[[] for _ in range(n)]
    for u,v in zip(src,dst):adj[u].append(v);rev[v].append(u)
    seen=set();order=[]
    for node in range(n):
        if node in seen:continue
        stack=[(node,False)]
        while stack:
            u,done=stack.pop()
            if done:order.append(u);continue
            if u in seen:continue
            seen.add(u);stack.append((u,True));stack.extend((v,False) for v in adj[u] if v not in seen)
    seen=set();components=[]
    for node in reversed(order):
        if node in seen:continue
        part=[];stack=[node];seen.add(node)
        while stack:
            u=stack.pop();part.append(u)
            for v in rev[u]:
                if v not in seen:seen.add(v);stack.append(v)
        components.append(sorted(part))
    components.sort(key=lambda x:(-len(x),x[0]));return components

def normalized(s):return re.sub(r'\s+',' ',s.strip().upper().replace('.',''))
def distance(a,b):
    lat1,lon1=a;lat2,lon2=b
    return math.hypot((lat1-lat2)*111195,(lon1-lon2)*111195*math.cos(math.radians((lat1+lat2)/2)))

def topology(geometry):
    endpoints=defaultdict(list)
    for sid,row in geometry.items():
        for side,field in [('start','from_street'),('end','to_street')]:
            key=tuple(sorted([normalized(row['street']),normalized(row[field])]))
            point=(float(row[side+'_latitude']),float(row[side+'_longitude']))
            endpoints[key].append((sid,side,point))
    nodes=[];assign={};splits=0
    for key,items in sorted(endpoints.items()):
        groups=[]
        for item in sorted(items):
            group=next((g for g in groups if all(distance(item[2],x[2])<=100 for x in g)),None)
            if group is None:group=[];groups.append(group)
            group.append(item)
        splits+=max(0,len(groups)-1)
        for group in groups:
            idx=len(nodes);coords=np.mean([x[2] for x in group],axis=0).tolist();nodes.append(dict(id=idx,names=key,lat=coords[0],lon=coords[1],endpoints=len(group)))
            for sid,side,point in group:assign[sid,side]=idx
    return nodes,assign,splits

def main():
    p=argparse.ArgumentParser();p.add_argument('--raw',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    root=Path(a.raw);out=Path(a.out)
    if out.exists():raise FileExistsError(out)
    begin=datetime.date(2019,1,1);dates=[begin+datetime.timedelta(days=i) for i in range(45)]
    records={r['date']:r for r in map(json.loads,(root/'download_manifest.jsonl').read_text(encoding='utf-8-sig').splitlines())}
    if set(records)!=set(map(str,dates)):raise ValueError('all 45 daily count checks must finish before preparation')
    variants=defaultdict(dict);source_files=[]
    for date in dates:
        rec=records[str(date)];assert rec['status'] in ('COMPLETE','NO_DATA')
        for item in rec['files']:
            path=root/item['name'];assert digest(path)==item['sha256'],'download checksum changed';source_files.append(item)
            if item['kind']=='geometry':
                with path.open(encoding='utf-8-sig',newline='') as f:
                    for row in csv.DictReader(f):
                        sid=int(row['segment_id']);key=tuple(row[k] for k in sorted(row) if k!='segment_id');variants[sid][key]=row
    unstable=sorted(s for s,v in variants.items() if len(v)>1);geometry={s:next(iter(v.values())) for s,v in variants.items() if len(v)==1}
    invalid=[]
    for sid,row in geometry.items():
        nums=[float(row[k]) for k in ['length','start_latitude','start_longitude','end_latitude','end_longitude']]
        if not all(map(math.isfinite,nums)) or nums[0]<=0 or not all(row[k].strip() for k in ['street','from_street','to_street']):invalid.append(sid)
    for sid in invalid:del geometry[sid]
    sids=sorted(geometry);columns={sid:e for e,sid in enumerate(sids)};speed=np.full((45,95,len(sids)),np.nan,dtype=np.float32)
    duplicates=conflicts=unknown=0;daily=[]
    for day,date in enumerate(dates):
        series=defaultdict(dict);rows=0
        for item in records[str(date)]['files']:
            if item['kind']!='speeds':continue
            with (root/item['name']).open(encoding='utf-8-sig',newline='') as f:
                for row in csv.DictReader(f):
                    rows+=1;sid=int(row['segment_id'])
                    if sid not in geometry:unknown+=1;continue
                    stamp=datetime.datetime.fromisoformat(row['time']);assert stamp.date()==date
                    t=(stamp-datetime.datetime.combine(date,datetime.time())).total_seconds();v=float(row['speed'])
                    if t in series[sid]:
                        duplicates+=1
                        if series[sid][t]!=v:series[sid][t]=float('nan');conflicts+=1
                    else:series[sid][t]=v
        assert rows==records[str(date)]['downloaded_rows']
        for sid,items in series.items():
            ordered=sorted(items);speed[day,:,columns[sid]]=aggregate_reports(ordered,[items[t] for t in ordered])[1:]
        daily.append(dict(date=str(date),source_rows=rows,complete_edge_quarters=int(np.isfinite(speed[day]).sum())))
        print('aggregated',date,rows,daily[-1]['complete_edge_quarters'],flush=True)
    valid=np.isfinite(speed);edge_keep=(valid.sum(axis=(0,1))>=40)&(valid.any(axis=1).sum(axis=0)>=5)
    retained={sid:geometry[sid] for sid,k in zip(sids,edge_keep) if k};nodes,assign,splits=topology(retained)
    chosen=sorted(retained);src=[assign[s,'start'] for s in chosen];dst=[assign[s,'end'] for s in chosen]
    components=largest_scc(len(nodes),src,dst)
    if not components:raise ValueError('No eligible road graph')
    keep=set(components[0]);chosen=[s for s in chosen if assign[s,'start'] in keep and assign[s,'end'] in keep and assign[s,'start']!=assign[s,'end']]
    index={v:i for i,v in enumerate(sorted(keep))};n=len(index);m=len(chosen)
    out.mkdir(parents=True);audit=dict(source_segments=len(variants),unstable_segments=unstable,invalid_geometry=invalid,
        coverage_eligible_edges=int(edge_keep.sum()),candidate_vertices=len(nodes),scc_sizes=[len(x) for x in components],n=n,m=m,
        same_name_coordinate_splits=splits,duplicate_records=duplicates,conflicting_records=conflicts,unknown_or_excluded_records=unknown,daily=daily)
    (out/'topology_audit.json').write_text(json.dumps(audit,indent=2))
    if n<100 or m<100:
        (out/'UNSUITABLE.txt').write_text('Largest connected directed graph has fewer than 100 vertices or edges. No traffic experiment was run.');print(json.dumps(audit));return
    matrix=speed[:,:, [columns[s] for s in chosen]].reshape(45*95,m);valid=np.isfinite(matrix)
    threshold=max(50,math.ceil(.1*m));low=(valid.sum(axis=1)<threshold);matrix[low]=np.nan;valid=np.isfinite(matrix)
    ww,ee=np.nonzero(valid);rng=np.random.default_rng(20261105);blocks=rng.random((45,m))<.05;held=blocks[ww//95,ee].astype(np.uint8)
    train=valid.copy();train[ww[held==1],ee[held==1]]=False;ff=np.array([np.quantile(matrix[train[:,e],e],.95) if train[:,e].any() else np.nan for e in range(m)],dtype='<f4')
    fallback=int((~np.isfinite(ff)).sum());ff[~np.isfinite(ff)]=np.nanmedian(ff)
    length=np.array([float(geometry[s]['length'])*1609.344 for s in chosen],dtype='<f4')
    with (out/'arcs.bin').open('wb') as f:
        np.array([index[assign[s,'start']] for s in chosen],dtype='<i4').tofile(f);np.array([index[assign[s,'end']] for s in chosen],dtype='<i4').tofile(f)
        length.tofile(f);ff.tofile(f);np.zeros(m,dtype='<i4').tofile(f)
    with (out/'coords.bin').open('wb') as f:
        np.array([nodes[v]['lon'] for v in sorted(keep)],dtype='<f4').tofile(f);np.array([nodes[v]['lat'] for v in sorted(keep)],dtype='<f4').tofile(f)
    np.array([d.weekday() for d in dates],dtype='<i4').tofile(out/'days.bin');np.arange(1,96,dtype='<i4').tofile(out/'slots.bin')
    obs=np.empty(len(ee),dtype=[('e','<i4'),('w','<i4'),('kph','<f4')]);obs['e']=ee;obs['w']=ww;obs['kph']=matrix[ww,ee];obs.tofile(out/'obs.bin')
    held.tofile(out/'holdout.bin');blocks.astype(np.uint8).tofile(out/'holdout_edge_day.bin')
    (out/'meta.txt').write_text(f'{n} {m} 45 96 900 0 95 {len(ee)} 1\n')
    with (out/'edge_mapping.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['edge','segment_id','u','v','length_m','estimated_freeflow_kph'])
        for e,s in enumerate(chosen):writer.writerow([e,s,index[assign[s,'start']],index[assign[s,'end']],length[e],ff[e]])
    with (out/'coverage.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['world','date','slot','valid_edges','retained'])
        for w in range(len(matrix)):writer.writerow([w,str(dates[w//95]),w%95+1,int(valid[w].sum()),int(not low[w])])
    info=dict(city='Chicago',source='https://data.cityofchicago.org/resource/sxs8-h27x',license='SEE_TERMS_OF_USE',n=n,m=m,days=45,dates=list(map(str,dates)),N=45*95,
        slot_seconds=900,slots_per_day=96,observations=len(ee),density=float(valid.mean()),heldout=int(held.sum()),heldout_fraction=float(held.mean()),
        holdout_unit='edge-day block',holdout_seed=20261105,freeflow='95th percentile of training aggregate speeds',freeflow_fallback_edges=fallback,
        reference='published CTA-based congestion estimates; not raw GPS trips',missing_policy='negative/zero invalid immediately; max 900-second report hold; full 900-second positive coverage required',
        aggregation='overlap-weighted arithmetic pace; first quarter-hour excluded to prevent source-report leakage across edge-day holdout blocks',
        empty_worlds=int(low.sum()),snapshot_min_edges=threshold,geometry_policy='stable segment geometry; named intersections only; complete-link 100 m maximum coordinate separation',
        no_data_dates=[r['date'] for r in records.values() if r['status']=='NO_DATA'],source_files=source_files,
        classes=['unknown'],one_snapshot_per_day_slot=True,upsampling=False,min_positive_kph=float(matrix[valid].min()),max_positive_kph=float(matrix[valid].max()))
    (out/'preparation.json').write_text(json.dumps(info,indent=2));print(json.dumps({k:v for k,v in info.items() if k not in ('source_files','dates')}),flush=True)
if __name__=='__main__':main()
