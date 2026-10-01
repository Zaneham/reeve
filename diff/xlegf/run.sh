#!/usr/bin/env bash
# Reeve differential for lib/xlegf.ml. Rebuilds both sides from scratch, sweeps
# the same inputs through each, and compares the output.
#
#   ./run.sh
#
# Run it from Git Bash on Windows: the OCaml side builds with the Windows opam
# and the Fortran side builds with gfortran inside WSL.
#
# Expect 2042 differences on Windows and none on Linux, because the Windows libm
# zeroes the low bits of cos near pi/2 and the Legendre series cancels onto it.
set -eu

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
work=$here/work
slatec=/c/dev/numerical/slatec-modern/src/original/src
opam=/c/Users/GGPC/AppData/Local/Microsoft/WinGet/Packages/OCaml.opam_Microsoft.Winget.Source_8wekyb3d8bbwe/opam.exe
w() { echo "$1" | sed 's,^/\([a-zA-Z]\)/,/mnt/\1/,'; }

rm -rf "$work"
mkdir -p "$work"
cp "$here/sweep.ml" "$work/"

echo "== building the library and the OCaml driver =="
(cd "$root" && "$opam" exec -- dune build)
objs=$root/_build/default/lib/.reeve.objs
(cd "$work" && "$opam" exec -- ocamlopt -I "$objs/byte" -I "$objs/native" \
  "$root/_build/default/lib/reeve.cmxa" -o sweep.exe sweep.ml >/dev/null)

echo "== generating the inputs and the OCaml answers =="
(cd "$work" && ./sweep.exe)

echo "== building the Fortran driver =="
src=""
for f in dxset dxadd dxadj dxred dxcon dxc210 dxlegf dxnrmp dxpmu dxpmup \
         dxpnrm dxpqnu dxpsi dxqmu dxqnu; do
  src="$src $(w "$slatec")/$f.f"
done
wwork=$(w "$work")
where=$(w "$here")
wsl -e bash -c "cd '$wwork' && gfortran -O0 -std=legacy -fdefault-integer-8 -w \
  '$where/mach.f90' '$where/xm.f' $src '$where/drv.f90' -o drv"

echo "== running the Fortran side =="
run_fortran() {
  local guard=$1 iguard=$2 dest=$3
  : > "$work/$dest"
  for f in "$work"/in.*; do
    wsl -e bash -c "cd '$wwork' && ./drv $guard $iguard < '$(basename "$f")' \
      >> '$dest' 2>/dev/null"
  done
}
run_fortran 0000000000000000 0 f.out
run_fortran 3FF0000000000000 7 f2.out

echo "== guard fill independence =="
if cmp -s "$work/f.out" "$work/f2.out"; then
  echo "the printed window does not depend on the guard fill"
else
  echo "THE PRINTED WINDOW DEPENDS ON THE GUARD FILL"
  cmp "$work/f.out" "$work/f2.out" || true
fi

echo "== cases and lines =="
cat "$work"/in.* | grep -cv '^[0-9]' | sed 's/^/cases:        /'
wc -l < "$work/ml.out" | sed 's/^/ml.out lines: /'
wc -l < "$work/f.out" | sed 's/^/f.out lines:  /'

echo "== comparison =="
if cmp "$work/f.out" "$work/ml.out"; then
  echo "IDENTICAL"
else
  echo "DIFFERENT"
  cp "$here/classify.py" "$work/"
  (cd "$work" && python classify.py) || true
  exit 1
fi
