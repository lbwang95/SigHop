#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
for r in 2 4 6 8; do
  bin/gen grid --k 10 --r "$r" --days 14 --slots 8 --oneway 0.3 --seed "$r" --out "$T/g.bin" >/dev/null
  for th in 1 8; do
    echo "ffref test r=$r build_threads=$th"
    bin/htq selftest --data "$T/g.bin" --tmp "$T/x.idx" --queries 80 --threads "$th" --pointsonly-tests 1 --ffref-tests 1
  done
done
echo "ALL FFREF TESTS PASSED"
