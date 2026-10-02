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

## Using it

    opam pin add reeve git+https://github.com/Zaneham/reeve

then `(libraries reeve)` in your dune file.

    open Reeve

    let info = Lapack.dgesv n 1 a 0 n ipiv 0 b 0 n
    let info = Lapack.dpotrf Blasmat.Upper n a 0 n
    let info = Lapack.dgetrf m n a 0 lda ipiv 0

    Blasmat.dgemm No_trans Trans m n k 1.0 a 0 lda b 0 ldb 0.0 c 0 ldc
    Blas.daxpy n alpha x 0 1 y 0 1
    let d = Blas.ddot n x 0 1 y 0 1

    let nz = Bessel.dbesj x 0.0 3 y
    let g  = Specfun.dgamln x
    let r  = Specfun.drf x y z
    let q, ierr, err = Quad.dgaus8 f a b 1e-12
    Sort.dsort x carry n Sort.Increasing_carry

Arguments are the Fortran's in the Fortran's order, with a base offset after
every array. `info` is returned, pivots are zero based, and character options
are variants.

Names are the Fortran's, lowercase. Each `.mli` has a line per routine.

## Build

    opam exec -- dune build
    opam exec -- dune test

## Licence

Apache 2.0

The original SLATEC library was released into the public domain by the Federal Government of the United States of America.
