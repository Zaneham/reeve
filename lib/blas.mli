(* Reeve, Copyright 2026 Zane Hambly.
   BLAS level 1, ported from the reference Fortran. A vector is a
   [float array] plus a base offset [xo] and an increment [incx], so the [i]
   th element lives at [xo + i * incx] when the increment is positive, and a
   negative increment walks the vector backwards from its far end as it does
   in the Fortran. *)

val daxpy : int -> float -> float array -> int -> int -> float array -> int
  -> int -> unit
(** [daxpy n da x xo incx y yo incy] adds [da] times the [n] element vector
    [x] from [xo] with stride [incx] to [y] from [yo] with stride [incy],
    doing nothing when [n] is not positive or [da] is zero. *)

val dscal : int -> float -> float array -> int -> int -> unit
(** [dscal n da x xo incx] scales the [n] elements of [x] from [xo] with
    stride [incx] by [da], doing nothing when [n] is not positive, [incx] is
    not positive or [da] is one. *)

val dcopy : int -> float array -> int -> int -> float array -> int -> int
  -> unit
(** [dcopy n x xo incx y yo incy] copies the [n] element vector [x] from [xo]
    with stride [incx] into [y] from [yo] with stride [incy]. *)

val dswap : int -> float array -> int -> int -> float array -> int -> int
  -> unit
(** [dswap n x xo incx y yo incy] exchanges the [n] elements of [x] from [xo]
    with stride [incx] against those of [y] from [yo] with stride [incy]. *)

val drot : int -> float array -> int -> int -> float array -> int -> int
  -> float -> float -> unit
(** [drot n x xo incx y yo incy c s] applies the plane rotation with cosine
    [c] and sine [s] to the [n] pairs taken from [x] from [xo] with stride
    [incx] and [y] from [yo] with stride [incy], overwriting both. *)

val drotg : float -> float -> float * float * float * float
(** [drotg da db] constructs the plane rotation that takes the vector
    [(da, db)] onto the first axis, returning [(r, z, c, s)] in that order:
    the rotated length [r], the single value [z] from which [c] and [s] can be
    reconstructed, and the cosine [c] and sine [s] themselves. *)

val ddot : int -> float array -> int -> int -> float array -> int -> int
  -> float
(** [ddot n x xo incx y yo incy] is the inner product of the [n] element
    vector [x] from [xo] with stride [incx] and [y] from [yo] with stride
    [incy], and 0 when [n] is not positive. The summation order is the
    Fortran's, including its path unrolled by five for unit strides. *)

val dasum : int -> float array -> int -> int -> float
(** [dasum n x xo incx] is the sum of the absolute values of the [n] elements
    of [x] from [xo] with stride [incx], and 0 when [n] is not positive or
    [incx] is not positive. The summation order is the Fortran's, including
    its path unrolled by six for a unit stride. *)

val dnrm2 : int -> float array -> int -> int -> float
(** [dnrm2 n x xo incx] is the Euclidean norm of the [n] element vector [x]
    from [xo] with stride [incx], computed by the scaled sum of squares so
    that it neither overflows nor underflows unnecessarily, and 0 when [n] is
    not positive. *)

val idamax : int -> float array -> int -> int -> int
(** [idamax n x xo incx] returns the index, counted from zero along the
    increment, of the first element of largest absolute value among the [n]
    elements of [x] starting at [xo] with stride [incx], and [-1] when [n] is
    below one or [incx] is not positive, which is where the Fortran returns
    its zero. *)
