#!/usr/bin/env python3
"""Check model and query files against data/manifest.json."""
import hashlib,json,struct
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def main():
    records=json.loads((ROOT/'data/manifest.json').read_text())
    for row in records:
        p=ROOT/row['path'];content=p.read_bytes()
        assert len(content)==row['bytes'] and hashlib.sha256(content).hexdigest()==row['sha256'],p
        if 'header' in row:
            assert content[:8]==b'HTQDATA1',p
            n,m,r,N,days,slots,coords,seconds=struct.unpack('<7qd',content[8:72])
            assert len(content)==72+8*m+4*m*r+4*N*r+12*N+(8*n if coords else 0),p
            assert list(row['header'].values())==[n,m,r,N,days,slots,coords,seconds],p
        else:assert len(content.splitlines())==row['queries'],p
    print('BUNDLED INPUTS VERIFIED:',len(records),'files')
if __name__=='__main__':main()
