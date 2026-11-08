"""Reader for the published HTQDATA1 inputs and history-prefix checks."""
import hashlib,struct
from pathlib import Path
import numpy as np
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
def load(path):
    with Path(path).open('rb') as f:
        assert f.read(8)==b'HTQDATA1'
        hdr=list(struct.unpack('<7q',f.read(56)));n,m,r,N,days,slots,coords=hdr
        seconds=struct.unpack('<d',f.read(8))[0]
        read=lambda kind,count:np.fromfile(f,dtype=kind,count=count)
        d=dict(hdr=hdr,seconds=seconds,fr=read('<i4',m),to=read('<i4',m),
               u=read('<f4',m*r).reshape(m,r),z=read('<f4',N*r).reshape(N,r),
               day=read('<i4',N),slot=read('<i4',N),dow=read('<i4',N),xy=f.read(n*8 if coords else 0))
        assert not f.read(1)
    assert seconds==900 and slots==96
    assert np.all(d['z'][:,0]==1)
    for k in ['u','z']:assert np.isfinite(d[k]).all() and (d[k]>=0).all()
    assert len(np.unique(d['day'].astype('i8')*96+d['slot']))==N
    return d
