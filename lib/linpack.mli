(* Reeve, Copyright 2026 Zane Hambly.
   LINPACK, ported from the SLATEC double precision routines. {!Raw} holds the
   Fortran shapes, where a matrix is a [float array] in column major order
   with an explicit leading dimension, so the element (i, j) zero based lives
   at [i + j * lda], and pivot vectors hold zero based row indices. Under that
   is the OCaml surface, which takes a {!Mat.t}, allocates its own pivots and
   workspace, and raises instead of handing back [info]. Band matrices are
   still a {!Mat.t}, since the band layout is a dense array of its own shape,
   but the bandwidths have to come in as arguments because the matrix does not
   carry them. *)

module Raw : sig
  (** The Fortran shapes, argument for argument. These are what the
      differential in [diff/] compares against the reference Fortran. Use the
      functions below unless you specifically want [lda]. *)

  val dgefa : float array -> int -> int -> int array -> int
  (** [dgefa a lda n ipvt] factors the order [n] matrix [a] by Gaussian
      elimination with partial pivoting, overwriting [a] with the multipliers
      below the diagonal and [u] on and above it, and filling [ipvt] with the
      [n] pivot rows. Returns [info], which is 0 normally and [k] when
      [u (k,k)] is zero, counting [k] from one as the Fortran does. *)

  val dgesl :
    float array -> int -> int -> int array -> float array -> int -> unit
  (** [dgesl a lda n ipvt b job] solves the order [n] system in place in [b],
      using the factors and pivots that {!dgefa} or {!dgeco} left in [a] and
      [ipvt]. [job] of 0 solves [a * x = b], any other value solves
      [transpose a * x = b]. *)

  val dgeco :
    float array -> int -> int -> int array -> float array -> float
  (** [dgeco a lda n ipvt z] factors the order [n] matrix [a] the way
      {!dgefa} does, filling [a] and [ipvt] the same way, and additionally
      estimates the condition of [a]. Returns [rcond], the reciprocal
      condition estimate, and leaves in the length [n] work vector [z] an
      approximate null vector when [rcond] is small. *)

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
      place in [b], using the factors and pivots that {!dgbfa} left in [abd]
      and [ipvt]. [job] of 0 solves [a * x = b], any other value solves
      [transpose a * x = b]. *)

  val dpofa : float array -> int -> int -> int
  (** [dpofa a lda n] factors the order [n] symmetric positive definite
      matrix [a] as [transpose r * r], reading and overwriting only the upper
      triangle. Returns [info], which is 0 on success and [k] when the leading
      order [k] minor is not positive definite, counting [k] from one as the
      Fortran does. *)

  val dposl : float array -> int -> int -> float array -> unit
  (** [dposl a lda n b] solves the order [n] system in place in [b], using the
      Cholesky factor that {!dpofa} left in [a]. *)

  val dpbfa : float array -> int -> int -> int -> int
  (** [dpbfa abd lda n m] factors the order [n] symmetric positive definite
      band matrix held in [abd] with [m] superdiagonals as
      [transpose r * r], overwriting [abd]. [abd] holds the diagonal in row
      [m] zero based and needs [lda] at least [m + 1]. Returns [info] with the
      same meaning as in {!dpofa}. *)

  val dpbsl : float array -> int -> int -> int -> float array -> unit
  (** [dpbsl abd lda n m b] solves the order [n] band system in place in [b],
      using the Cholesky factor that {!dpbfa} left in [abd]. *)
end

(* ---- The OCaml surface ---- *)

exception Singular of int
(** A pivot came out exactly zero at this step, counting from one. *)

exception Not_positive_definite of int
(** The leading minor of this order is not positive definite. *)

val solve : Mat.t -> Mat.t -> unit
(** [solve a b] solves [a x = b] for every column of [b], overwriting [b]
    with [x] and [a] with its factors. Raises {!Singular}. *)

val lu : Mat.t -> int array
(** [lu a] factors [a] in place and returns the pivots, which are zero based
    row indices, one less than the Fortran's. Raises {!Singular}. *)

val lu_solve : ?trans:bool -> Mat.t -> int array -> Mat.t -> unit
(** [lu_solve a ipvt b] solves using the factors {!lu} or {!lu_rcond} left
    behind, overwriting every column of [b]. [~trans:true] solves with the
    transpose. *)

val lu_rcond : Mat.t -> int array * float
(** [lu_rcond a] factors [a] the way {!lu} does and also estimates its
    condition, returning the pivots and the reciprocal estimate. A singular
    [a] gives you 0.0 rather than an exception, since the estimate is the
    point. *)

val band_solve : ml:int -> mu:int -> Mat.t -> Mat.t -> unit
(** [band_solve ~ml ~mu abd b] solves the band system for every column of
    [b]. [abd] is the band layout, with [a i j] at row [ml + mu + i - j] and
    at least [2 * ml + mu + 1] rows, the extra [ml] of them scratch for the
    pivoting. Raises {!Singular}. *)

val band_lu : ml:int -> mu:int -> Mat.t -> int array
(** [band_lu ~ml ~mu abd] factors the band matrix in place and returns the
    pivots. Same layout as {!band_solve}, and it raises {!Singular} the same
    way. *)

val band_lu_solve :
  ?trans:bool -> ml:int -> mu:int -> Mat.t -> int array -> Mat.t -> unit
(** [band_lu_solve ~ml ~mu abd ipvt b] solves using the factors {!band_lu}
    left behind, overwriting every column of [b]. [~trans:true] solves with
    the transpose. *)

val cholesky : Mat.t -> unit
(** [cholesky a] factors the symmetric positive definite [a] as
    [transpose r * r] in place, reading and writing only the upper triangle.
    Raises {!Not_positive_definite}. *)

val cholesky_solve : Mat.t -> Mat.t -> unit
(** [cholesky_solve a b] solves using the factor {!cholesky} left behind,
    overwriting every column of [b]. *)

val solve_spd : Mat.t -> Mat.t -> unit
(** [solve_spd a b] factors and solves in one go for symmetric positive
    definite [a]. Raises {!Not_positive_definite}. *)

val band_cholesky : mu:int -> Mat.t -> unit
(** [band_cholesky ~mu abd] factors the symmetric positive definite band
    matrix in place. [abd] holds [a i j] for [i <= j] at row [mu + i - j] and
    needs [mu + 1] rows. Raises {!Not_positive_definite}. *)

val band_cholesky_solve : mu:int -> Mat.t -> Mat.t -> unit
(** [band_cholesky_solve ~mu abd b] solves using the factor
    {!band_cholesky} left behind, overwriting every column of [b]. *)

val band_solve_spd : mu:int -> Mat.t -> Mat.t -> unit
(** [band_solve_spd ~mu abd b] factors and solves in one go for a symmetric
    positive definite band matrix. Raises {!Not_positive_definite}. *)
