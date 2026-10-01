(* Reeve, Copyright 2026 Zane Hambly.
   LAPACK's dense linear system routines, the LU and Cholesky paths, ported
   from the reference Fortran. A matrix is a [float array] plus a base offset
   [ao] and a leading dimension [lda], column major, so the element (i, j)
   zero based lives at [ao + i + j * lda]. A pivot vector is an [int array]
   plus a base offset and holds zero based row indices, one less than the
   value the Fortran stores. Every routine returns the Fortran's [info]: zero
   on success and a positive index when the factor is numerically singular or
   not positive definite. A bad argument, where the Fortran calls [xerbla],
   raises [Invalid_argument] naming the routine and the Fortran argument
   position instead. *)

val dgetf2 : int -> int -> float array -> int -> int -> int array -> int
  -> int
(** [dgetf2 m n a ao lda ipiv ipo] factors the [m] by [n] matrix [a] as
    [p * l * u] by unblocked Gaussian elimination with partial pivoting,
    overwriting [a] and filling [min m n] zero based pivot rows into [ipiv]
    from [ipo]. Returns [info], which is [k] when [u (k,k)] is exactly zero,
    counting [k] from one. *)

val dgetrf2 : int -> int -> float array -> int -> int -> int array -> int
  -> int
(** [dgetrf2 m n a ao lda ipiv ipo] is the recursive right looking version of
    {!dgetf2}, with the same arguments, the same output and the same [info].
    This is what {!dgetrf} calls. *)

val dgetrf : int -> int -> float array -> int -> int -> int array -> int
  -> int
(** [dgetrf m n a ao lda ipiv ipo] factors the [m] by [n] matrix [a] as
    [p * l * u] by the blocked algorithm, calling {!dgetrf2} on each panel,
    with the same output and the same [info] as {!dgetf2}. *)

val dgetrs : Blasmat.trans -> int -> int -> float array -> int -> int
  -> int array -> int -> float array -> int -> int -> int
(** [dgetrs trans n nrhs a ao lda ipiv ipo b bo ldb] solves the order [n]
    system with [nrhs] right hand sides in place in the [n] by [nrhs] matrix
    [b], using the factors and zero based pivots that {!dgetrf} left in [a]
    and [ipiv]. [trans] of [No_trans] solves [a * x = b], anything else
    solves [transpose a * x = b]. Returns [info], always zero here. *)

val dgesv : int -> int -> float array -> int -> int -> int array -> int
  -> float array -> int -> int -> int
(** [dgesv n nrhs a ao lda ipiv ipo b bo ldb] solves [a * x = b] for the
    order [n] matrix [a] and [nrhs] right hand sides, factoring [a] in place
    with {!dgetrf} and overwriting [b] with the solution. Returns [info],
    which is [k] when [u (k,k)] is exactly zero, in which case [b] is left
    alone. *)

val dpotf2 : Blasmat.uplo -> int -> float array -> int -> int -> int
(** [dpotf2 uplo n a ao lda] factors the order [n] symmetric positive
    definite matrix [a] as [transpose u * u] when [uplo] is [Blasmat.Upper]
    and [l * transpose l] when it is [Blasmat.Lower], by the unblocked
    algorithm, reading and overwriting only that triangle. Returns [info],
    which is [k] when the leading order [k] minor is not positive definite,
    counting [k] from one. *)

val dpotrf2 : Blasmat.uplo -> int -> float array -> int -> int -> int
(** [dpotrf2 uplo n a ao lda] is the recursive version of {!dpotf2}, with the
    same arguments, the same output and the same [info]. This is what
    {!dpotrf} calls. *)

val dpotrf : Blasmat.uplo -> int -> float array -> int -> int -> int
(** [dpotrf uplo n a ao lda] factors the order [n] symmetric positive
    definite matrix [a] by the blocked algorithm, calling {!dpotrf2} on each
    diagonal block, with the same output and the same [info] as
    {!dpotf2}. *)

val dpotrs : Blasmat.uplo -> int -> int -> float array -> int -> int
  -> float array -> int -> int -> int
(** [dpotrs uplo n nrhs a ao lda b bo ldb] solves the order [n] system with
    [nrhs] right hand sides in place in the [n] by [nrhs] matrix [b], using
    the Cholesky factor that {!dpotrf} left in the [uplo] triangle of [a].
    Returns [info], always zero here. *)

val dposv : Blasmat.uplo -> int -> int -> float array -> int -> int
  -> float array -> int -> int -> int
(** [dposv uplo n nrhs a ao lda b bo ldb] solves [a * x = b] for the order
    [n] symmetric positive definite matrix [a] and [nrhs] right hand sides,
    factoring [a] in place with {!dpotrf} and overwriting [b] with the
    solution. Returns [info], which is [k] when the leading order [k] minor
    is not positive definite, in which case [b] is left alone. *)
