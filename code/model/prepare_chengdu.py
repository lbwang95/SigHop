#!/usr/bin/env python3
"""Adapt the public Chengdu link/speed files, keeping their actual 2-minute time tags."""
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from graph_utils import largest_scc

p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--out',required=True)
p.add_argument('--seed',type=int,default=20261104);a=p.parse_args()
root=Path(a.root);out=Path(a.out)
if out.exists():raise FileExistsError(out)
catalog=json.loads((root/'catalog.json').read_text())
with (root/'link.csv').open(encoding='utf-8-sig',newline='') as f:links=list(csv.DictReader(f))
ids=sorted({int(r['Node_Start']) for r in links}|{int(r['Node_End']) for r in links})
index={x:i for i,x in enumerate(ids)}
src=np.array([index[int(r['Node_Start'])] for r in links]);dst=np.array([index[int(r['Node_End'])] for r in links])
keep=largest_scc(len(ids),src,dst);kept_edges=keep[src]&keep[dst]
chosen=[r for r,k in zip(links,kept_edges) if k];nodes=[x for x,k in zip(ids,keep) if k]
index={x:i for i,x in enumerate(nodes)};n=len(nodes);m=len(chosen)
edge_map={int(r['Link']):e for e,r in enumerate(chosen)}
length=np.array([float(r['Length']) for r in chosen],dtype='<f4')
if not (np.isfinite(length).all() and (length>0).all()):raise RuntimeError('invalid published length')
coords={}
for r in links:
    coords.setdefault(int(r['Node_Start']),(float(r['Longitude_Start']),float(r['Latitude_Start'])))
    coords.setdefault(int(r['Node_End']),(float(r['Longitude_End']),float(r['Latitude_End'])))
slot_tags=np.array([s for lo,hi in [(90,150),(240,300),(360,420),(510,570),(630,690)] for s in range(lo,hi)],dtype='<i4')
slot_index={int(s):i for i,s in enumerate(slot_tags)}
speed_items=[f for f in catalog['files'] if f['name'].startswith('speed')]
dates=sorted({datetime.date(2015,int(re.search(r'\[(\d+)\]',f['name'])[1])//100,int(re.search(r'\[(\d+)\]',f['name'])[1])%100) for f in speed_items})
day_index={d:i for i,d in enumerate(dates)}
speed=np.full((m,len(dates)*300),np.nan,dtype=np.float32)
excluded=[];source_records=[]
for item in speed_items:
    meta=json.loads((root/f'{item["id"]}.json').read_text())
    file=root/meta['packed_file']
    with file.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
    if digest!=meta['packed_sha256']:raise RuntimeError('packed file checksum differs')
    datecode,half=map(int,re.findall(r'\[(\d+)\]',item['name']))
    date=datetime.date(2015,datecode//100,datecode%100)
    with np.load(file,allow_pickle=False) as z:
        link=z['link'];slot=z['slot'];values=z['speed']
    expected=set(map(int,slot_tags[:150] if half==0 else slot_tags[150:]))
    actual=set(map(int,np.unique(slot)))
    if not actual.issubset(expected):
        excluded.append(dict(id=item['id'],name=item['name'],reason='periods disagree with filename half',source_md5=item['computed_md5']))
        continue
    keys=slot.astype(np.int64)*100000+link
    if np.unique(keys).size!=len(keys):raise RuntimeError('duplicate link/period within one source file')
    ee=np.fromiter((edge_map.get(int(x),-1) for x in link),dtype=np.int32,count=len(link))
    ww=np.fromiter((day_index[date]*300+slot_index[int(x)] for x in slot),dtype=np.int32,count=len(slot))
    select=ee>=0
    speed[ee[select],ww[select]]=values[select]
    source_records.append(dict(id=item['id'],name=item['name'],md5=item['computed_md5'],rows=len(link)))
valid=np.isfinite(speed)&(speed>0)
ww,ee=np.nonzero(valid.T)
rng=np.random.default_rng(a.seed);held=(rng.random(len(ee))<.05).astype(np.uint8)
training=valid.copy();training[ee[held==1],ww[held==1]]=False
ff=np.full(m,np.nan,dtype=np.float32)
for e in range(m):
    if training[e].any():ff[e]=np.quantile(speed[e,training[e]],.95)
fallback=int(np.isnan(ff).sum());ff[~np.isfinite(ff)]=np.nanmedian(ff)
if not (np.isfinite(ff).all() and (ff>0).all()):raise RuntimeError('freeflow unavailable')
out.mkdir(parents=True,exist_ok=False)
with (out/'arcs.bin').open('wb') as f:
    np.array([index[int(r['Node_Start'])] for r in chosen],dtype='<i4').tofile(f)
    np.array([index[int(r['Node_End'])] for r in chosen],dtype='<i4').tofile(f)
    length.tofile(f);ff.astype('<f4').tofile(f);np.zeros(m,dtype='<i4').tofile(f)
with (out/'coords.bin').open('wb') as f:
    np.array([coords[x][0] for x in nodes],dtype='<f4').tofile(f)
    np.array([coords[x][1] for x in nodes],dtype='<f4').tofile(f)
np.array([d.weekday() for d in dates],dtype='<i4').tofile(out/'days.bin')
slot_tags.tofile(out/'slots.bin');held.tofile(out/'holdout.bin')
obs=np.empty(len(ee),dtype=[('e','<i4'),('w','<i4'),('kph','<f4')])
obs['e']=ee;obs['w']=ww;obs['kph']=speed[ee,ww];obs.tofile(out/'obs.bin')
(out/'meta.txt').write_text(f'{n} {m} {len(dates)} 720 120 0 300 {len(obs)} 1\n')
with (out/'edge_mapping.csv').open('w',newline='') as f:
    writer=csv.writer(f);writer.writerow(['edge','source_link','u','v','length_m','estimated_freeflow_kph'])
    for e,r in enumerate(chosen):writer.writerow([e,r['Link'],r['Node_Start'],r['Node_End'],r['Length'],float(ff[e])])
report=dict(city='Chengdu',source='https://doi.org/10.6084/m9.figshare.7140209.v4',license=catalog['license'],
            source_nodes=len(ids),source_edges=len(links),n=n,m=m,dropped_vertices=len(ids)-n,
            days=len(dates),dates=list(map(str,dates)),N=len(dates)*300,slot_seconds=120,
            time_windows=['03:00-05:00','08:00-10:00','12:00-14:00','17:00-19:00','21:00-23:00'],
            observations=len(obs),density=float(valid.mean()),heldout=int(held.sum()),
            empty_worlds=int((valid.sum(axis=0)==0).sum()),excluded_source_files=excluded,
            freeflow='95th percentile of training positive published speeds per edge',freeflow_fallback_edges=fallback,
            classes=['unknown (not supplied by source)'],reference='published processed speed data; raw GPS observation mask unavailable',
            min_positive_kph=float(speed[valid].min()),max_positive_kph=float(speed[valid].max()),seed=a.seed,
            source_files=source_records,shortest_path_scope='largest strongly connected component')
(out/'preparation.json').write_text(json.dumps(report,indent=2))
print(json.dumps({k:v for k,v in report.items() if k not in ('dates','source_files')}),flush=True)
