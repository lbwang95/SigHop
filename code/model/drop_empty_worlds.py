#!/usr/bin/env python3
"""Remove wholly unobserved worlds after fitting, without changing any fitted value on retained worlds."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import struct
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--out',required=True);a=p.parse_args()
source=Path(a.input);out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
with (source/'raw.bin').open('rb') as f:
    assert f.read(8)==b'HTQRAW01';N,nnz=struct.unpack('<2q',f.read(16))
    offsets=np.fromfile(f,dtype='<i8',count=N+1)
    kept=np.flatnonzero(offsets[1:]>offsets[:-1])
    with (out/'raw.bin').open('wb') as g:
        g.write(b'HTQRAW01');g.write(struct.pack('<2q',len(kept),nnz))
        offsets[np.r_[kept,N]].astype('<i8').tofile(g)
        shutil.copyfileobj(f,g,length=1024*1024)
first_tags=None
for file in sorted(source.glob('dataset_r*.bin')):
    with file.open('rb') as f:
        assert f.read(8)==b'HTQDATA1';header=list(struct.unpack('<7q',f.read(56)));seconds=f.read(8)
        n,m,r,oldN,days,slots,coords=header
        assert oldN==N
        arcs=f.read(m*8);loadings=f.read(m*r*4)
        z=np.fromfile(f,dtype='<f4',count=N*r).reshape(N,r)
        tags=[np.fromfile(f,dtype='<i4',count=N) for _ in range(3)]
        xy=f.read(n*8 if coords else 0)
        assert not f.read(1)
    header[3]=len(kept)
    with (out/file.name).open('wb') as f:
        f.write(b'HTQDATA1');f.write(struct.pack('<7q',*header));f.write(seconds)
        f.write(arcs);f.write(loadings);z[kept].tofile(f)
        for values in tags:values[kept].tofile(f)
        f.write(xy)
    if first_tags is None:first_tags=tags
with (out/'world_mapping.csv').open('w',newline='') as f:
    writer=csv.writer(f);writer.writerow(['new_world','source_world','day_index','actual_slot','weekday'])
    for i,w in enumerate(kept):writer.writerow([i,int(w),*[int(x[w]) for x in first_tags]])
with (out/'removed_worlds.csv').open('w',newline='') as f:
    writer=csv.writer(f);writer.writerow(['source_world','day_index','actual_slot','weekday'])
    for w in np.flatnonzero(offsets[1:]==offsets[:-1]):writer.writerow([int(w),*[int(x[w]) for x in first_tags]])
for name in ('nmf_report.csv','edge_mapping.csv'):
    shutil.copyfile(source/name,out/name)
info=json.loads((source/'preparation.json').read_text())
info.update(N=len(kept),source_grid_worlds=N,empty_worlds_removed=N-len(kept),empty_worlds=0,
            density=nnz/(m*len(kept)),fit_input_directory=str(source),
            world_filter_rule='omit only worlds with zero published positive speed entries; fitted U and retained Z unchanged')
(out/'preparation.json').write_text(json.dumps(info,indent=2))
print(json.dumps({'source_worlds':N,'retained_worlds':len(kept),'removed':N-len(kept),'speed_entries':nnz}))
