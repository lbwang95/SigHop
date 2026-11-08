#!/usr/bin/env bash
# Builds the binaries and checks every method against one Dijkstra search per snapshot on small
# grids, with 1 and 8 build threads. The last line must be ALL TESTS PASSED.
set -euo pipefail
cd "$(dirname "$0")/.."
make -s
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
./bin/gen grid --k 10 --r 2 --days 14 --slots 8 --seed 2 --out "$T/g2.bin" >/dev/null
./bin/gen grid --k 12 --r 4 --days 14 --slots 8 --seed 4 --out "$T/g4.bin" >/dev/null
./bin/gen grid --k 18 --r 6 --days 30 --slots 12 --oneway 0.3 --seed 11 --out "$T/src6.bin" >/dev/null
./bin/gen grid --k 18 --r 6 --days 60 --z-from "$T/src6.bin" --seed 12 --out "$T/boot6.bin" >/dev/null
./bin/gen grid --k 16 --r 8 --days 20 --slots 10 --oneway 0.3 --seed 21 --out "$T/g8.bin" >/dev/null
for f in g2 g4 boot6 g8; do
  for th in 1 8; do
    echo "== $f threads=$th"
    ./bin/htq selftest --data "$T/$f.bin" --tmp "$T/x.idx" --queries 80 --threads $th | tail -2 | head -1
    ./bin/htq selftest --data "$T/$f.bin" --tmp "$T/x.idx" --queries 80 --threads $th >/dev/null
  done
done
echo "ALL TESTS PASSED"
