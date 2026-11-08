#!/usr/bin/env python3
"""Download the frozen DIMACS inputs outside Git, checking every SHA-256."""
import argparse,hashlib,json,shutil,tarfile,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(8<<20),b''):h.update(chunk)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--archive',type=Path,help='use an already downloaded archive');a=ap.parse_args()
    meta=json.loads((ROOT/'data/dimacs_manifest.json').read_text())
    archive=a.archive or ROOT/'data'/meta['asset']
    if not archive.exists():
        if a.archive:raise FileNotFoundError(archive)
        temp=archive.with_suffix(archive.suffix+'.part')
        with urllib.request.urlopen(meta['url'],timeout=60) as response,temp.open('wb') as output:shutil.copyfileobj(response,output,8<<20)
        if sha(temp)!=meta['sha256']:raise ValueError('Archive SHA-256 mismatch')
        temp.replace(archive)
    if sha(archive)!=meta['sha256']:raise ValueError('Archive SHA-256 mismatch')
    expected={r['path']:r for r in meta['files']}
    with tarfile.open(archive,'r:gz') as tar:
        members=tar.getmembers()
        if len(members)!=len(expected) or {m.name for m in members}!=set(expected):raise ValueError('Unexpected archive members')
        for member in members:
            if not member.isfile():raise ValueError('Archive member is not a regular file')
            path=ROOT/'data'/member.name;record=expected[member.name]
            if path.exists():
                if sha(path)!=record['sha256']:raise FileExistsError('Different existing input: '+str(path))
                continue
            path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix(path.suffix+'.part')
            with tar.extractfile(member) as source,temp.open('wb') as out:shutil.copyfileobj(source,out,8<<20)
            if temp.stat().st_size!=record['bytes'] or sha(temp)!=record['sha256']:raise ValueError('Input checksum mismatch: '+member.name)
            temp.replace(path)
    print('DIMACS INPUTS VERIFIED (8 files)')
if __name__=='__main__':main()
