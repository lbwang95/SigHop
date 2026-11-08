#!/usr/bin/env python3
"""Download public CC0 CSVs, verify original MD5, and losslessly package the float32 input representation."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import time
import urllib.request
import numpy as np
import pyarrow.csv as pc

p=argparse.ArgumentParser();p.add_argument('--catalog',default='chengdu_catalog.json');p.add_argument('--out',default='Chengdu')
p.add_argument('--limit',type=int,default=0);a=p.parse_args()
catalog=json.loads(Path(a.catalog).read_text())
assert catalog['license']['name']=='CC0'
root=Path(a.out);raw=root/'raw';packed=root/'packed'
raw.mkdir(parents=True,exist_ok=True);packed.mkdir(exist_ok=True)
shutil.copyfile(a.catalog,packed/'catalog.json')
for item in catalog['files'][:a.limit or None]:
    name=item['name'];source=raw/name;meta=packed/(str(item['id'])+'.json')
    if meta.exists():continue
    start=time.monotonic()
    if not source.exists():
        for attempt in range(3):
            try:
                with urllib.request.urlopen(item['download_url'],timeout=60) as r,source.with_suffix('.part').open('wb') as f:
                    while block:=r.read(1024*1024):f.write(block)
                source.with_suffix('.part').replace(source);break
            except Exception:
                if attempt==2:raise
                time.sleep(3)
    with source.open('rb') as f: md5=hashlib.file_digest(f,'md5').hexdigest()
    if md5!=item['computed_md5'] or source.stat().st_size!=item['size']:raise RuntimeError(f'checksum differs: {name}')
    record=dict(item);record['source_md5_verified']=True
    if name=='link.csv':
        shutil.copyfile(source,packed/name)
        record['packed_file']=name
    else:
        assert re.fullmatch(r'speed\[\d+\]_\[[01]\]\.csv',name)
        table=pc.read_csv(source,read_options=pc.ReadOptions(use_threads=False))
        period=table['Period'].combine_chunks().dictionary_encode()
        slots=[]
        for label in period.dictionary.to_pylist():
            begin=label.split('-')[0];hh,mm=map(int,begin.split(':'))
            if mm%2:raise RuntimeError('unexpected period')
            slots.append(hh*30+mm//2)
        slot=np.array(slots,dtype='<i2')[period.indices.to_numpy()]
        link=table['Link'].to_numpy().astype('<i4')
        speed=table['Speed'].to_numpy().astype('<f4')
        target=packed/(str(item['id'])+'.npz')
        np.savez_compressed(target,link=link,slot=slot,speed=speed)
        with target.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
        record.update(packed_file=target.name,packed_sha256=digest,packed_bytes=target.stat().st_size,
                      rows=len(speed),slots=sorted(set(slots)),finite_positive=int((np.isfinite(speed)&(speed>0)).sum()),
                      representation='float32 speed, int32 link id, int16 actual 2-minute slot; original CSV retained')
    record['seconds']=time.monotonic()-start
    meta.write_text(json.dumps(record,indent=2))
    print(json.dumps({k:record[k] for k in ('id','name','size','seconds')}),flush=True)
print('PACKING COMPLETE',flush=True)
