#!/usr/bin/env python3
"""Validate the X2 prefix contract as a scheduler-generated preparation job."""
import csv, sys
from pathlib import Path
import numpy as np
from data_io import load,sha

out=Path('results/provenance/x2_prefix.csv');out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w',newline='') as f:
    writer=csv.writer(f);writer.writerow(['graph','short_N','long_N','checks','mismatches','short_sha256','long_sha256'])
    for g in sys.argv[1:]:
        short=Path(f'../data/generated/{g}_r4_d90.bin');long=Path(f'../data/generated/{g}_r4_d365.bin')
        a,b=load(short),load(long);n=a['hdr'][3]
        comparisons=[a['hdr'][i]==b['hdr'][i] for i in (0,1,2,5,6)]
        comparisons += [a['seconds']==b['seconds'],a['xy']==b['xy']]
        comparisons += [np.array_equal(a[k],b[k]) for k in ('fr','to','u')]
        comparisons += [np.array_equal(a[k],b[k][:n]) for k in ('z','day','slot','dow')]
        mismatches=sum(not v for v in comparisons)
        writer.writerow([g,n,b['hdr'][3],len(comparisons),mismatches,sha(short),sha(long)]);f.flush()
        print(f'{g} prefix checks={len(comparisons)} mismatches={mismatches}',flush=True)
        if mismatches:raise RuntimeError('X2 prefix contract failed; STOP')
