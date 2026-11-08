#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
for r in 2 4 6 8; do
  bin/gen grid --k 10 --r "$r" --days 14 --slots 8 --oneway 0.3 --seed "$r" --out "$T/g.bin" > /dev/null
  for th in 1 8; do
    echo "== pointsonly r=$r threads=$th (in-memory and reloaded, all query variants)"
    bin/htq selftest --data "$T/g.bin" --tmp "$T/test.idx" --queries 80 --threads "$th" --pointsonly-tests 1
  done
done
echo "ALL POINTSONLY TESTS PASSED"
