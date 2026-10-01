# Reeve

A numerical library for OCaml. It's a port of SLATEC, with a few LAPACK
routines, and it's written in OCaml throughout, so the only thing it depends on
is the compiler.

A reeve was the manor official who checked the tallies were true in case you are curious.
This is also the sister project to [zblas](https://github.com/Zaneham/zblas) which converts a few of these 
routines to IBM high level assembler.

The project these routines have been ported from is [SLATEC](https://github.com/Zaneham/SLATEC).
This is my modernised version, you can find the original at
<https://www.netlib.org/slatec/>.

## Modules

`Blas` and `Blasmat` are the reference BLAS, levels 1 to 3. `Lapack` is the LU
and Cholesky paths, with `Laux` under it. `Linpack` is SLATEC's LINPACK.
`Specfun`, `Bessel`, `Gamma`, `Expint` and `Xlegf` are the special functions.
`Quad` is quadrature, `Interp` is polynomial interpolation, `Sort` is the
sorting routines.

## Checking

`diff/` runs the same inputs through the reference Fortran and through Reeve
and compares IEEE patterns. BLAS and LAPACK agree over 15,736 cases and
1,453,154 lines of output. Sorting agrees over 409,800 lines, the exponential
integrals over 2,897 values, quadrature over 56 cases.

Coefficient tables are pulled out of the Fortran by script and compared back by
the pattern, not read by eye.

Level 4 of my modernised SLATEC scheme asks what hostile gfortran flags do to a
routine. OCaml has no equivalent, so Reeve does not claim it. However as this was made to benchmark
some of my stuff in the OCaml compiler I may write some numbers down

## Calling conventions

A matrix is a `float array`, a base offset and an `lda`, column major, so
`(i, j)` is at `off + i + j * lda`. A vector is a `float array`, a base offset
and an increment. Fortran character options are variants. Pivot and permutation
vectors are zero based.

## Build

    opam exec -- dune build
    opam exec -- dune test

## Licence

Apache 2.0

The original SLATEC library was released into the public domain by the Federal Government of the United States of America.
