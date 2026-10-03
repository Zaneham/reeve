#!/bin/sh
# Reeve, Copyright 2026 Zane Hambly.
# Times the bench under two opam switches, alternating between them so that
# thermal and frequency drift cancels instead of landing on whichever ran
# second. Reports the best of all rounds per switch.
#
#   sh bench/compare.sh <switch-a> <switch-b> [rounds]

set -e
A=${1:?switch a}
B=${2:?switch b}
R=${3:-5}
OPAM=/c/Users/GGPC/AppData/Local/Microsoft/WinGet/Packages/OCaml.opam_Microsoft.Winget.Source_8wekyb3d8bbwe/opam.exe
HERE=$(cd "$(dirname "$0")/.." && pwd)
cd "$HERE"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

for s in "$A" "$B"; do
  "$OPAM" exec --switch "$s" -- dune build --build-dir "_build.$s" bench/bench.exe >/dev/null 2>&1
done

i=1
while [ "$i" -le "$R" ]; do
  for s in "$A" "$B"; do
    "$OPAM" exec --switch "$s" -- dune exec --build-dir "_build.$s" bench/bench.exe \
      2>/dev/null | awk 'NR>1 && NF>=5 {print $1" "$3" "$5}' >> "$OUT/$s"
  done
  printf 'round %s of %s\n' "$i" "$R" >&2
  i=$((i + 1))
done

awk -v a="$A" -v b="$B" '
FILENAME ~ a"$" { if (!(($1) in pa) || $2+0 < pa[$1]) pa[$1] = $2+0; ca[$1] = $3; next }
                { if (!(($1) in pb) || $2+0 < pb[$1]) pb[$1] = $2+0; cb[$1] = $3; ord[++k] = $1 }
END {
  printf "%-13s %12s %12s %8s  %s\n", "bench", a, b, "change", "checksum"
  for (i = 1; i <= k; i++) {
    n = ord[i]
    printf "%-13s %12.1f %12.1f %+7.2f%%  %s\n", n, pa[n], pb[n],
      (pb[n] - pa[n]) / pa[n] * 100.0, (ca[n] == cb[n] ? "match" : "DIFFER")
  }
}' "$OUT/$A" "$OUT/$B"
