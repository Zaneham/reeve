(* Reeve, Copyright 2026 Zane Hambly.
   LINPACK, ported from the SLATEC double precision routines. A matrix is a
   [float array] in column major order with an explicit leading dimension, so
   the element (i, j) zero based lives at [i + j * lda]. Pivot vectors hold
   zero based row indices. *)

val dgefa : float array -> int -> int -> int array -> int
(** [dgefa a lda n ipvt] factors the order [n] matrix [a] by Gaussian
    elimination with partial pivoting, overwriting [a] with the multipliers
    below the diagonal and [u] on and above it, and filling [ipvt] with the
    [n] pivot rows. Returns [info], which is 0 normally and [k] when [u (k,k)]
    is zero, counting [k] from one as the Fortran does. *)

val dgesl : float array -> int -> int -> int array -> float array -> int -> unit
(** [dgesl a lda n ipvt b job] solves the order [n] system in place in [b],
    using the factors and pivots that {!dgefa} or {!dgeco} left in [a] and
    [ipvt]. [job] of 0 solves [a * x = b], any other value solves
    [transpose a * x = b]. *)

val dgeco : float array -> int -> int -> int array -> float array -> float
(** [dgeco a lda n ipvt z] factors the order [n] matrix [a] the way {!dgefa}
    does, filling [a] and [ipvt] the same way, and additionally estimates the
    condition of [a]. Returns [rcond], the reciprocal condition estimate, and
    leaves in the length [n] work vector [z] an approximate null vector when
    [rcond] is small. *)

val dgbfa : float array -> int -> int -> int -> int -> int array -> int
(** [dgbfa abd lda n ml mu ipvt] factors the order [n] band matrix held in
    [abd] with [ml] subdiagonals and [mu] superdiagonals, overwriting [abd]
    and filling [ipvt] with the [n] pivot rows. [abd] holds the diagonal in
    row [ml + mu] zero based and needs [lda] at least [2 * ml + mu + 1].
    Returns [info] with the same meaning as in {!dgefa}. *)

val dgbsl :
  float array -> int -> int -> int -> int -> int array -> float array -> int
  -> unit
(** [dgbsl abd lda n ml mu ipvt b job] solves the order [n] band system in
    place in [b], using the factors and pivots that {!dgbfa} left in [abd] and
    [ipvt]. [job] of 0 solves [a * x = b], any other value solves
    [transpose a * x = b]. *)

val dpofa : float array -> int -> int -> int
(** [dpofa a lda n] factors the order [n] symmetric positive definite matrix
    [a] as [transpose r * r], reading and overwriting only the upper triangle.
    Returns [info], which is 0 on success and [k] when the leading order [k]
    minor is not positive definite, counting [k] from one as the Fortran
    does. *)

val dposl : float array -> int -> int -> float array -> unit
(** [dposl a lda n b] solves the order [n] system in place in [b], using the
    Cholesky factor that {!dpofa} left in [a]. *)

val dpbfa : float array -> int -> int -> int -> int
(** [dpbfa abd lda n m] factors the order [n] symmetric positive definite band
    matrix held in [abd] with [m] superdiagonals as [transpose r * r],
    overwriting [abd]. [abd] holds the diagonal in row [m] zero based and
    needs [lda] at least [m + 1]. Returns [info] with the same meaning as in
    {!dpofa}. *)

val dpbsl : float array -> int -> int -> int -> float array -> unit
(** [dpbsl abd lda n m b] solves the order [n] band system in place in [b],
    using the Cholesky factor that {!dpbfa} left in [abd]. *)
