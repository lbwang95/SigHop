#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
out="${1:-results/example}"
if [[ -e "$out" ]]; then echo "Choose a fresh output directory: $out" >&2; exit 1; fi
mkdir -p "$out"
python3 scripts/verify_data.py
make
bin/htq build --data ../data/models/CD_r4.bin --name CD_r4 --prune points \
  --threads 8 --index "$out/CD_r4.idx" --csv "$out/build.csv"
for method in sighop cchbatch cch dijkstra ttp; do
  extra=()
  if [[ "$method" == sighop ]]; then extra+=(--index "$out/CD_r4.idx"); fi
  bin/htq query --data ../data/models/CD_r4.bin --name CD_r4 --method "$method" \
    --qdir ../data/queries/CD --sets D8 --per 20 --rounds 1 --warmup 0 \
    --phi 'cd1730:dow=0-4;days=15-45;slots=70-71' \
    --timeout 60 --skip-after 5 "${extra[@]}" --csv "$out/query.csv"
done
python3 - "$out/query.csv" <<'PY'
import csv,sys
rows=list(csv.DictReader(open(sys.argv[1])))
assert len(rows)==5
assert all(int(r['completed'])==20 and int(r['mismatches'])==0 for r in rows)
print('EXAMPLE PASSED: all five methods, 20 queries each, mismatches=0')
PY
