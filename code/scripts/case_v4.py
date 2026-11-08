#!/usr/bin/env python3
"""Chengdu case analysis using deterministic shortest-path searches."""
import csv,hashlib,heapq,json,math,struct,time
from collections import defaultdict,deque
from pathlib import Path
import numpy as np

R=Path('results/real');P=Path('results/provenance')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(path,header,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:w=csv.writer(f);w.writerow(header.split(','));w.writerows(rows)

def load(path):
    with open(path,'rb') as f:
        assert f.read(8)==b'HTQDATA1'
        n,m,r,N,days,slots,coords=struct.unpack('<7q',f.read(56));seconds=struct.unpack('<d',f.read(8))[0]
        fr=np.fromfile(f,'<i4',m);to=np.fromfile(f,'<i4',m)
        u=np.fromfile(f,'<f4',m*r).reshape(m,r).astype('f8');z=np.fromfile(f,'<f4',N*r).reshape(N,r).astype('f8')
        day,slot,dow=[np.fromfile(f,'<i4',N) for _ in range(3)]
        xy=[np.fromfile(f,'<f4',n) for _ in range(2)] if coords else [None,None]
        assert not f.read(1)
    assert (n,m,r,N,days,slots,seconds)==(1902,5943,4,1780,45,96,900)
    assert np.isfinite(u).all() and np.isfinite(z).all() and (u>=0).all() and (z>=0).all()
    assert np.all(z[:,0]==1) and len(np.unique(day.astype('i8')*96+slot))==N
    return n,m,r,N,fr,to,u,z,day,slot,dow,xy

def main():
    R.mkdir(parents=True,exist_ok=True);P.mkdir(parents=True,exist_ok=True)
    assert not (R/'v4_case_profile.csv').exists(),'Output already exists: results/real/v4_case_profile.csv'
    data=Path('../data/models/CD_r4.bin');old=Path('../results/v2/real/case_cd.csv')
    inputs=[data,old]+[Path(f'../data/queries/CD/D{i}.txt') for i in [4,6,8,10]]
    hashes={str(p):sha(p) for p in inputs};begin=time.time()
    n,m,r,N,fr,to,u,z,day,slot,dow,xy=load(data)
    # Same arc order as GraphCsr: original edge id ascending within each outgoing list.
    adj=[[] for _ in range(n)]
    for e,(a,b) in enumerate(zip(fr,to)):adj[int(a)].append((int(b),e))
    def dijkstra(source,weights):
        dist=[math.inf]*n;parent=[-1]*n;dist[source]=0.0;heap=[(0.0,source)]
        while heap:
            du,x=heapq.heappop(heap)
            if du>dist[x]:continue
            for y,e in adj[x]:
                nd=du+weights[e]
                if nd<dist[y]:
                    dist[y]=nd;parent[y]=e;heapq.heappush(heap,(nd,y))
        return np.asarray(dist),np.asarray(parent,dtype='i4')
    def path(parent,s,t):
        arcs=[];v=t
        while v!=s:
            e=int(parent[v]);assert e>=0 and len(arcs)<n,(s,t,v)
            assert int(to[e])==v
            arcs.append(e);v=int(fr[e])
        arcs.reverse();return tuple(arcs)
    weights={};selected={};mean_weights={};world_lists={}
    available=sorted(set(map(int,slot)))
    assert len(available)==40
    for tau in available:
        ids=np.flatnonzero((day>=15)&(day<45)&(dow>=0)&(dow<=4)&(slot==tau))
        assert len(ids)==22 and len(set(map(int,day[ids])))==22
        selected[tau]=ids
        # Preserve legacy dotd component accumulation order. Inputs originate as float32.
        w=np.zeros((len(ids),m),dtype='f8')
        for j in range(r):w+=z[ids,j,None]*u[None,:,j]
        assert (w>0).all() and np.isfinite(w).all()
        weights[tau]=w;world_lists[tau]=[row.tolist() for row in w]
        mu=z[ids].mean(axis=0);wm=np.zeros(m,dtype='f8')
        for j in range(r):wm+=u[:,j]*mu[j]
        mean_weights[tau]=wm.tolist()
    def compute(s,tau):
        results=[dijkstra(s,w) for w in world_lists[tau]]
        ds=np.stack([x[0] for x in results]);parents=np.stack([x[1] for x in results])
        mean,parent=dijkstra(s,mean_weights[tau]);sorted_dist=np.sort(ds,axis=0)
        return ds,parents,mean,parent,sorted_dist
    def qr(alpha,k):return max(0,min(k-1,math.ceil(alpha*k-1e-9)-1))
    def metrics(s,t,tau,cache):
        ds,parents,mean,parent,ordered=cache;k=len(selected[tau]);assert np.isfinite(ds[:,t]).all()
        q=[float(ordered[qr(a,k),t]) for a in [.5,.8,.95]]
        fixed=path(parent,s,t);fv=np.zeros(k,dtype='f8')
        for e in fixed:fv+=weights[tau][:,e]
        assert np.all(fv+1e-8*np.maximum(1,fv)>=ds[:,t])
        fv.sort();fq=[float(fv[qr(a,k)]) for a in [.5,.8,.95]]
        routes={path(parents[i],s,t) for i in range(k)}
        return q,float(mean[t]),int(np.sum(ds[:,t]>mean[t])),fq,len(routes)
    start=lambda tau:f'{tau//4:02d}:{tau%4*15:02d}'
    anchor={};profile=[];validation=[]
    reference={x['start']:x for x in csv.DictReader(old.open())}
    for tau in available:
        cache=compute(1291,tau);anchor[tau]=cache
        q,eta,late,fq,routes=metrics(1291,1127,tau,cache)
        row=[start(tau),22,*[x/60 for x in q],eta/60,late,*[x/60 for x in fq],routes];profile.append(row)
        expected=reference[start(tau)];errors=[abs(q[i]/60-float(expected[f'p{percent}_min']))/max(1e-12,abs(float(expected[f'p{percent}_min']))) for i,percent in enumerate([50,80,95])]
        mism=int(expected['days'])!=22 or routes!=int(expected['routes']) or max(errors)>1e-4
        validation.append([start(tau),22,*errors,routes,int(expected['routes']),int(mism)])
    write(P/'v4_case_validation.csv','start,days,p50_relative_error,p80_relative_error,p95_relative_error,routes,baseline_routes,mismatches',validation)
    mismatches=sum(x[-1] for x in validation);print(f'case profile mismatches={mismatches}/40',flush=True)
    if mismatches:raise RuntimeError('CD profile differs from the reference results')
    write(R/'v4_case_profile.csv','start,days,p50_min,p80_min,p95_min,mean_eta_min,late_mean,fixed_p50_min,fixed_p80_min,fixed_p95_min,routes',profile)
    queries=defaultdict(list);unique_counts={}
    for setno in [4,6,8,10]:
        pairs=[tuple(map(int,line.split())) for line in Path(f'../data/queries/CD/D{setno}.txt').read_text().splitlines() if line.strip()]
        assert len(pairs)==1000 and all(0<=s<n and 0<=t<n and s!=t for s,t in pairs)
        unique_counts[f'D{setno}']=len(set(pairs))
        if setno==8:assert pairs[0]==(1291,1127)
        for qi,(s,t) in enumerate(pairs):queries[s].append((setno,qi,t))
    pairs_out=[]
    for number,s in enumerate(sorted(queries),1):
        for tau in [12,32,48,70,84]:
            cache=anchor[tau] if s==1291 else compute(s,tau)
            for setno,qi,t in queries[s]:
                q,eta,late,fq,routes=metrics(s,t,tau,cache)
                pairs_out.append((setno,qi,tau,[f'D{setno}',s,t,start(tau),22,*q,eta,late,fq[0],fq[2],routes]))
        if number%50==0:print('case sources',number,'/',len(queries),flush=True)
    pairs_out.sort(key=lambda x:x[:3]);assert len(pairs_out)==20000
    write(R/'v4_case_pairs.csv','set,s,t,start,days,p50,p80,p95,mean_eta,late_mean,fixed_p50,fixed_p95,routes',[x[3] for x in pairs_out])
    reachable={1291};queue=deque([1291])
    while queue:
        for t,e in adj[queue.popleft()]:
            if t not in reachable:reachable.add(t);queue.append(t)
    targets=sorted(reachable-{1291});assert len(targets)==1901
    maps=[]
    for tau in [12,70]:
        ds,parents,mean,parent,ordered=anchor[tau]
        for t in targets:maps.append([t,float(xy[0][t]) if xy[0] is not None else 'NA',float(xy[1][t]) if xy[1] is not None else 'NA',start(tau),float(ordered[qr(.5,22),t]),float(ordered[qr(.95,22),t]),float(mean[t])])
    write(R/'v4_case_map.csv','t,x,y,start,p50,p95,mean_eta',maps)
    query_file=Path('../data/queries/CD/MAP.txt')
    expected_map=''.join(f'1291 {t}\n' for t in targets)
    if query_file.exists():assert query_file.read_text()==expected_map, 'Existing MAP queries differ'
    else:query_file.write_text(expected_map)
    (P/'v4_MAP.txt').write_bytes(query_file.read_bytes())
    assert all(sha(p)==h for p,h in hashes.items())
    protocol=dict(inputs=hashes,query_unique_counts=unique_counts,distinct_sources=len(queries),source=1291,targets=len(targets),profile_rows=40,pair_rows=20000,map_rows=len(maps),days=22,day_range=[15,45],slots=available,
        ties='float64 sequential component sums; min-heap (distance,vertex-id); outgoing arcs in original edge-id order; predecessor changes only for strictly smaller distance (same as legacy pathQuery)',
        mean_eta='shortest path under mean factor vector, not mean of daily optimal distances',late_mean='count of days daily optimal T_i exceeds mean_eta; lower bound on lateness if the mean-condition route is held fixed',
        map_timing='1901 independent point-to-point queries',map_query_sha256=sha(query_file),wall_seconds=time.time()-begin)
    (P/'v4_case_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    print('CASE_ANALYSIS_COMPLETE',json.dumps({k:protocol[k] for k in ['targets','profile_rows','pair_rows','map_rows','distinct_sources','wall_seconds']}),flush=True)

if __name__=='__main__':main()
