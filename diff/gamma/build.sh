#!/bin/sh
# Builds both sides of the gamma differential from scratch and compares them.
#
#   ./build.sh            native mingw gfortran, which shares libm with OCaml
#   ./build.sh wsl        gfortran inside WSL instead
#
# The reference Fortran is compiled out of the slatec-modern tree rather than
# vendored here.  -finit-integer=0 pins the unassigned i that d9gmic, d9gmit,
# d9lgic and d9lgit test their convergence guard against, so the guard never
# fires, which is what the port does inside the convergent region.
# -finit-real=zero pins dgamit's algap1 and sgngam on the branch where the
# Fortran leaves them unassigned, to the 0.0 and 1.0 the port passes.
set -e

here=$(cd "$(dirname "$0")" && pwd)
slatec=/c/dev/numerical/slatec-modern
reeve=$(cd "$here/../.." && pwd)
opam=/c/Users/GGPC/AppData/Local/Microsoft/WinGet/Packages/OCaml.opam_Microsoft.Winget.Source_8wekyb3d8bbwe/opam.exe
flags="-O2 -finit-integer=0 -finit-real=zero"
out="$here/build"

mkdir -p "$out"
rm -f "$out"/*.o "$out"/*.mod "$out"/*.cm* "$out"/*.exe "$out"/*.ml \
      "$out"/fdrv "$out"/fort.txt "$out"/ocaml.txt

if [ "$1" = wsl ]; then
  w=$(printf %s "$here" | sed 's|^/\([a-zA-Z]\)|/mnt/\1|')
  ws=$(printf %s "$slatec" | sed 's|^/\([a-zA-Z]\)|/mnt/\1|')
  wsl -e bash -c "set -e; cd '$w'; \
    gfortran $flags -J build -I '$ws/src/modern/special_functions' -c '$ws/src/modern/service/service.f90' -o build/service.o; \
    gfortran $flags -J build -I build -c cases.f90 -o build/cases.o; \
    gfortran $flags -J build -I build -c mach.f90 -o build/mach.o; \
    gfortran $flags -J build -I build -I '$ws/src/modern/special_functions' -c gmod.f90 -o build/gmod.o; \
    gfortran $flags -std=legacy -J build -c '$ws/src/original/src/dpsifn.f' -o build/dpsifn.o; \
    gfortran $flags -J build -I build -c fdrv.f90 -o build/fdrv.o; \
    gfortran $flags -o build/fdrv build/fdrv.o build/gmod.o build/cases.o build/mach.o build/dpsifn.o build/service.o; \
    ./build/fdrv > build/fort.txt"
else
  gf=/c/MinGW/bin/gfortran
  inc="$slatec/src/modern/special_functions"
  $gf $flags -J "$out" -I "$inc" -c "$slatec/src/modern/service/service.f90" -o "$out/service.o"
  $gf $flags -J "$out" -I "$out" -c "$here/cases.f90" -o "$out/cases.o"
  $gf $flags -J "$out" -I "$out" -c "$here/mach.f90" -o "$out/mach.o"
  $gf $flags -J "$out" -I "$out" -I "$inc" -c "$here/gmod.f90" -o "$out/gmod.o"
  $gf $flags -std=legacy -J "$out" -c "$slatec/src/original/src/dpsifn.f" -o "$out/dpsifn.o"
  $gf $flags -J "$out" -I "$out" -c "$here/fdrv.f90" -o "$out/fdrv.o"
  $gf $flags -o "$out/fdrv.exe" "$out/fdrv.o" "$out/gmod.o" "$out/cases.o" \
    "$out/mach.o" "$out/dpsifn.o" "$out/service.o"
  "$out/fdrv.exe" > "$out/fort.txt"
fi

"$opam" exec -- dune build --root "$reeve" @@default

# The OCaml native toolchain here is mingw and wants Windows paths, and
# ocamlopt leaves its .cmi and .cmx beside the source, so it compiles a copy.
objs=$(cygpath -w "$reeve/_build/default/lib/.reeve.objs")
cp "$here/odrv.ml" "$out/odrv.ml"
"$opam" exec -- ocamlopt \
  -I "$objs\\byte" -I "$objs\\native" \
  "$(cygpath -w "$reeve/_build/default/lib/reeve.cmxa")" \
  -o "$(cygpath -w "$out/odrv.exe")" "$(cygpath -w "$out/odrv.ml")" >/dev/null
"$out/odrv.exe" > "$out/ocaml.txt"

echo "fortran $(wc -l < "$out/fort.txt") lines, ocaml $(wc -l < "$out/ocaml.txt") lines"
if cmp "$out/fort.txt" "$out/ocaml.txt"; then
  echo "identical"
else
  echo "differ, see build/fort.txt and build/ocaml.txt"
fi
