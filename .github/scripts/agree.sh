#!/usr/bin/env bash
set -euo pipefail

name=${1:?usage: agree.sh NAME FORTRAN-OUTPUT REEVE-OUTPUT MINIMUM-LINES}
fort=${2:?usage: agree.sh NAME FORTRAN-OUTPUT REEVE-OUTPUT MINIMUM-LINES}
ml=${3:?usage: agree.sh NAME FORTRAN-OUTPUT REEVE-OUTPUT MINIMUM-LINES}
floor=${4:?usage: agree.sh NAME FORTRAN-OUTPUT REEVE-OUTPUT MINIMUM-LINES}

for f in "$fort" "$ml"; do
  if [ ! -f "$f" ]; then
    echo "$name: $f was never written, so nothing was compared"
    exit 1
  fi
done

nf=$(wc -l < "$fort")
nm=$(wc -l < "$ml")

if [ "$nf" -eq 0 ] || [ "$nm" -eq 0 ]; then
  echo "$name: a driver produced no output, so nothing was compared"
  echo "$name: reference Fortran $nf lines, Reeve $nm lines"
  exit 1
fi

if [ "$nf" -lt "$floor" ] || [ "$nm" -lt "$floor" ]; then
  echo "$name: the sweep is short of its $floor line floor, so it covered less than it should"
  echo "$name: reference Fortran $nf lines, Reeve $nm lines"
  exit 1
fi

if [ "$nf" -ne "$nm" ]; then
  echo "$name: the two sides printed different numbers of lines"
  echo "$name: reference Fortran $nf lines, Reeve $nm lines"
  { diff "$fort" "$ml" || true; } | head -40
  exit 1
fi

if ! cmp -s "$fort" "$ml"; then
  differing=$({ diff "$fort" "$ml" || true; } | grep -c '^<' || true)
  echo "$name: $differing of $nf lines disagree"
  { diff "$fort" "$ml" || true; } | head -40
  exit 1
fi

echo "$name: $nf lines, identical"
