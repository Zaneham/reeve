#!/bin/sh
# Reeve differential, interpolation. Rebuilds both sides from scratch, sweeps
# the same inputs through the reference Fortran and through Reeve, and compares
# the two bit pattern dumps. Run from Git Bash in diff/interp.
#
#   ./run.sh
#
# gfortran lives in WSL, the OCaml toolchain on the Windows side.

#
# SLATEC_SRC overrides where the reference Fortran lives, as a WSL path.
set -e

HERE=$(cd "$(dirname "$0")" && pwd)
REEVE=$HERE/../..
WHERE=/mnt/c${HERE#/c}
OPAM=/c/Users/GGPC/AppData/Local/Microsoft/WinGet/Packages/OCaml.opam_Microsoft.Winget.Source_8wekyb3d8bbwe/opam.exe

rm -f "$HERE"/*.o "$HERE"/*.mod "$HERE"/*.cmi "$HERE"/*.cmx "$HERE"/drv.exe \
  "$HERE"/odrv.exe "$HERE"/cases.hex "$HERE"/f.out "$HERE"/o.out "$HERE"/f.err

wsl -e bash -lc "cd $WHERE && S=${SLATEC_SRC:-/mnt/c/dev/numerical/slatec-modern/src} && set -e \
  && gfortran -c -std=f2018 -O2 -ffp-contract=off \$S/modern/service/service.f90 -I\$S/modern/service \
  && gfortran -c -std=f2018 -O2 -ffp-contract=off mach.f90 \
  && gfortran -c -std=f2018 -O2 -ffp-contract=off slat.f90 -I\$S/modern/interpolation \
  && gfortran -c -std=legacy -O2 -ffp-contract=off \$S/original/src/dplint.f \$S/original/src/dpolcf.f \$S/original/src/dpolvl.f \
  && gfortran -o drv.exe drv.f90 slat.o mach.o service.o dplint.o dpolcf.o dpolvl.o \
    -std=f2018 -O2 -ffp-contract=off -I. \
  && ./drv.exe > f.out 2> f.err"

cat "$HERE/f.err"

( cd "$REEVE" && "$OPAM" exec -- dune build )
"$OPAM" exec -- ocamlopt \
  -I "$REEVE/_build/default/lib/.reeve.objs/byte" \
  -I "$REEVE/_build/default/lib/.reeve.objs/native" \
  "$REEVE/_build/default/lib/reeve.cmxa" "$HERE/odrv.ml" -o "$HERE/odrv.exe"

( cd "$HERE" && ./odrv.exe > o.out )

wc -l "$HERE/f.out" "$HERE/o.out"
if cmp "$HERE/f.out" "$HERE/o.out"; then
  echo IDENTICAL
else
  echo DIFFER
  exit 1
fi
