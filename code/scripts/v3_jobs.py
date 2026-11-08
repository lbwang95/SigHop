"""V3-only additions to the existing make_jobs.py entry point."""
import csv,sys
from pathlib import Path

def main():
    stage=sys.argv[sys.argv.index('--v3')+1]
    tasks=[]
    def job(jid,kind,mem,threads,deps,cmd):
        tasks.append('\t'.join([jid,kind,str(round(mem)),str(threads),','.join(deps) or '-',cmd]))
    def build(g,prune):
        nm=g+'_r4_d90';data='../data/generated/'+nm+'.bin'
        out='idx/NY_r4_d90_pointsonly.idx' if prune=='pointsonly' else 'idx/tmp_v3_'+nm+'_points.idx'
        memory=80 if prune=='pointsonly' else {'NY':40,'BAY':60,'COL':150,'FLA':300}[g]
        timeout=7200 if prune=='pointsonly' else 14400
        job('v3_build_'+nm+'_'+prune,'build',memory,32,[data],
            f'bin/htq build --data {data} --name {nm} --prune {prune} --index {out} --timeout {timeout} --threads 32 --csv results/build_v3.csv --stats-csv results/prune_stats.csv')
    if stage=='gate':build('NY','points')
    elif stage=='build':
        for g in ['BAY','COL','FLA']:build(g,'points')
        build('NY','pointsonly')
    elif stage=='time':
        for g in ['NY','BAY','COL','FLA']:
            nm=g+'_r4_d90';data='../data/generated/'+nm+'.bin';memory=8 if g=='FLA' else 7
            for w in [7,14,30,60]:
                phi=f"'wk{w}:dow=0-4;slots=32-33;days={90-w}-90'"
                job(f'e2_{nm}_ttp_w{w}','time',memory,1,[data,f'../data/queries/{g}/D10.txt'],
                    f'bin/htq query --data {data} --name {nm} --method ttp --qdir ../data/queries/{g} --sets D10 --phi {phi} --per 20 --rounds 1 --warmup 0 --timeout 60 --skip-after 5 --csv results/e2_ttp.csv')
        for k in [1,3,5]:
            job(f'v3_x6_NY_ttpk{k}','time',7,1,['../data/generated/NY_r4_d90.bin','../data/queries/NY/D10.txt'],
                f"bin/htq query --data ../data/generated/NY_r4_d90.bin --name NY_r4_d90_ttpk{k} --method ttp --qdir ../data/queries/NY --sets D1,D2,D3,D4,D10 --phi 'wk90:dow=0-4;slots=32-33;days=0-90' --per 20 --rounds 1 --warmup 0 --timeout 60 --skip-after 5 --ttp-k {k} --csv results/x6.csv")
        builds=list(csv.DictReader(open('results/build_v3.csv'))) if Path('results/build_v3.csv').exists() else []
        only=next((r for r in builds if r['name']=='NY_r4_d90' and r['prune']=='pointsonly' and r['status']=='OK'),None)
        if only:
            mem=next(r for r in csv.DictReader(open('../results/v2/mem.csv')) if r['name']=='NY_r4_d90')
            reserve=float(only['label_mb'])/1000*1.3+65*float(mem['cch_world_mb'])/1000+6
            sets=','.join('D%d'%i for i in range(1,11))
            job('v3_e8_NY_pointsonly','time',reserve,1,['../data/generated/NY_r4_d90.bin','../data/queries/NY/D10.txt','idx/NY_r4_d90_pointsonly.idx'],
                f"bin/htq query --data ../data/generated/NY_r4_d90.bin --name NY_r4_d90 --method sighop --index idx/NY_r4_d90_pointsonly.idx --qdir ../data/queries/NY --sets {sets} --phi 'wk90:dow=0-4;slots=32-33;days=0-90' --per 1000 --rounds 3 --warmup 16 --csv results/e8_pointsonly.csv")
    else:raise SystemExit('Unknown V3 stage: '+stage)
    print('\n'.join(tasks))
