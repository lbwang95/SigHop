"""Jobs for parameter baselines, ablations, and the Chengdu case study."""
import csv,sys
from pathlib import Path
GRAPHS=['NY','BAY','COL','FLA'];SETS=','.join('D%d'%i for i in range(1,11))
ALPHAS={1:'0.5',2:'0.5,0.95',5:'0.1,0.25,0.5,0.75,0.9',10:'0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95'}
PY=sys.executable

def read(p):
    if not Path(p).exists():return []
    with Path(p).open() as f:
        rd=csv.DictReader(f);return [r for r in rd if r!=dict(zip(rd.fieldnames,rd.fieldnames))]

def main():
    stage=sys.argv[sys.argv.index('--v4')+1];jobs=[]
    builds={(r['name'],r['prune']):r for p in ['results/build.csv','results/v4_build.csv'] for r in read(p) if r['status']=='OK'}
    mems={r['name']:r for r in read('results/mem.csv')}
    def job(jid,kind,mem,threads,deps,command):jobs.append('\t'.join([jid,kind,str(round(mem)),str(threads),','.join(deps) or '-',command]))
    def phi(tag,days=90,slot=32,dow='0-4',lo=0):return f"'{tag}:dow={dow};slots={slot}-{slot+1};days={lo}-{days}'"
    def kdays(days,dow=(0,1,2,3,4)):return sum(d%7 in dow for d in range(days))
    def tjob(exp,g,r,d,method,sets='D10',condition=None,k=None,alphas=None,tag='',prune='points',per=None,real=False):
        nm='CD_r4' if real else f'{g}_r{r}_d{d}';data='../data/models/CD_r4.bin' if real else f'../data/generated/{nm}.bin'
        idx='idx/CD_r4_points.idx' if real else f'idx/{nm}_{prune}.idx'
        if exp=='e8':idx='idx/tmp_v4_FLA_r4_d90_box.idx'
        m=mems[nm];world=float(m['cch_world_mb'])/1000;ver=k*world
        deps=[data,f'../data/queries/{g}/{"MAP" if real else "D10"}.txt'];extra=''
        if method.startswith('sighop'):
            b=builds.get((nm,prune))
            if not b:return
            reserve=float(b['label_mb'])/1000*1.3+ver+6;extra+=' --index '+idx;deps.append(idx)
        elif method in ('cch','cchbatch'):
            if float(m['cch_all_mb'])>1e6:return
            reserve=ver*(2 if method=='cchbatch' else 1)+6
        elif method=='ttp':reserve=k*float(m['n'])*4/1e9+ver+6
        else:reserve=ver+6
        per=per if per is not None else 20 if method=='ttp' else 100 if method=='dijkstra' else 1000
        extra+=f' --per {per}'
        if method=='ttp':extra+=' --rounds 1 --warmup 0 --timeout 60 --skip-after 5'
        elif method=='dijkstra':extra+=' --rounds 1 --warmup 1'
        else:extra+=' --rounds 3 --warmup 16'
        if alphas:extra+=' --alphas '+alphas
        output='results/real/v4_case_time.csv' if real else f'results/v4_{exp}.csv'
        job(f'v4_{exp}_{nm}_{method}{tag}','time',reserve,1,deps,
            f'bin/htq query --data {data} --name {nm} --method {method} --qdir ../data/queries/{g} --sets {sets} --phi {condition}{extra} --csv {output}')
    if stage=='default-check':
        tjob('default_check','NY',4,90,'sighop',sets='D1,D10',condition=phi('wk90'),k=65)
    elif stage.startswith('build-'):
        number=int(stage.split('-')[1]);spec=[('NY','points',14400),('NY','pointsonly',14400),('FLA','pointsonly',14400),('BAY','pointsonly',14400),('COL','pointsonly',14400),('FLA','box',7200)]
        if 1<=number<=6:
            g,prune,timeout=spec[number-1];nm=f'{g}_r4_d90';memory={'NY':40,'BAY':60,'COL':150,'FLA':300}[g]*(1 if prune=='points' else 2)
            job(f'v4_build_{number}_{nm}_{prune}','build',memory,32,[f'../data/generated/{nm}.bin'],
                f'bin/htq build --data ../data/generated/{nm}.bin --name {nm} --prune {prune} --timeout {timeout} --threads 32 --stats-csv results/v4_prune_stats.csv --csv results/v4_build.csv --index idx/tmp_v4_{nm}_{prune}.idx')
        elif number==7:
            if not Path('idx/CD_r4_points.idx').exists():
                job('v4_build_CD_r4','build',8,8,['../data/models/CD_r4.bin'],
                    'bin/htq build --data ../data/models/CD_r4.bin --name CD_r4 --prune points --threads 8 --index idx/CD_r4_points.idx --csv results/v4_build.csv')
        else:raise SystemExit('invalid build stage')
    elif stage=='case':
        job('v4_case_analysis','gen',16,1,['../data/models/CD_r4.bin']+[f'../data/queries/CD/D{i}.txt' for i in [4,6,8,10]],'make bin/model_case && '+PY+' scripts/case_v4.py')
    elif stage=='time':
        for label,slot,dow in [('wk03',12,'0-4'),('wk12',48,'0-4'),('wk17',68,'0-4'),('wk21',84,'0-4'),('we08',32,'5-6'),('all08',32,'0-6')]:
            k=65 if dow=='0-4' else 25 if dow=='5-6' else 90
            for method in ['dijkstra','ttp']:tjob('e9','NY',4,90,method,condition=phi(label,slot=slot,dow=dow),k=k,tag='_'+label)
        for n,alphas in ALPHAS.items():
            for method in ['cch','dijkstra','ttp']:tjob('e4','NY',4,90,method,condition=phi('wk90'),k=65,alphas=alphas,tag='_a%d'%n)
        for r in [2,6,8]:tjob('e3','NY',r,90,'ttp',condition=phi('wk90'),k=65)
        for d in [7,30,180,365,730,1095]:tjob('e7','NY',4,d,'ttp',condition=phi('wk%d'%d,days=d),k=kdays(d))
        for g in GRAPHS:tjob('e5',g,4,90,'sighop-ffref',sets=SETS,condition=phi('wk90'),k=65)
        for g in ['BAY','COL']:
            for method in ['sighop-nohub','sighop-countsep']:tjob('e5',g,4,90,method,sets=SETS,condition=phi('wk90'),k=65)
        if builds.get(('FLA_r4_d90','box')):
            tjob('e8','FLA',4,90,'sighop',sets='D1,D10',condition=phi('wk90'),k=65,prune='box')
        pairs=sum(bool(s.strip()) for s in Path('../data/queries/CD/MAP.txt').read_text().splitlines())
        for method in ['sighop','cchbatch','cch','dijkstra','ttp']:
            tjob('case','CD',4,45,method,sets='MAP',condition=phi('cd_wk1730',days=45,slot=70,lo=15),k=22,per=pairs,real=True)
    else:raise SystemExit('unknown V4 stage')
    print('\n'.join(jobs))
