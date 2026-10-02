(* Reeve, Copyright 2026 Zane Hambly.
   BLAS level 1, ported from the reference Fortran. The Fortran shapes live in
   {!Raw}, where a vector is a [float array] plus a base offset [xo] and an
   increment [incx], so the [i] th element lives at [xo + i * incx] when the
   increment is positive, and a negative increment walks the vector backwards
   from its far end as it does in the Fortran. Below that is the surface an
   OCaml caller wants, which takes whole arrays and works out the length
   itself. *)

module Raw : sig
  (** The Fortran shapes, argument for argument. These are what the
      differential in [diff/] compares against the reference Fortran, and
      they're what {!Lapack} and {!Linpack} are written against. Use the
      functions below unless you want an interior section of an array. *)

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
  (** [dcopy n x xo incx y yo incy] copies the [n] element vector [x] from
      [xo] with stride [incx] into [y] from [yo] with stride [incy]. *)

  val dswap : int -> float array -> int -> int -> float array -> int -> int
    -> unit
  (** [dswap n x xo incx y yo incy] exchanges the [n] elements of [x] from
      [xo] with stride [incx] against those of [y] from [yo] with stride
      [incy]. *)

  val drot : int -> float array -> int -> int -> float array -> int -> int
    -> float -> float -> unit
  (** [drot n x xo incx y yo incy c s] applies the plane rotation with cosine
      [c] and sine [s] to the [n] pairs taken from [x] from [xo] with stride
      [incx] and [y] from [yo] with stride [incy], overwriting both. *)

  val drotg : float -> float -> float * float * float * float
  (** [drotg da db] constructs the plane rotation that takes the vector
      [(da, db)] onto the first axis, returning [(r, z, c, s)] in that order:
      the rotated length [r], the single value [z] from which [c] and [s] can
      be reconstructed, and the cosine [c] and sine [s] themselves. *)

  val ddot : int -> float array -> int -> int -> float array -> int -> int
    -> float
  (** [ddot n x xo incx y yo incy] is the inner product of the [n] element
      vector [x] from [xo] with stride [incx] and [y] from [yo] with stride
      [incy], and 0 when [n] is not positive. The summation order is the
      Fortran's, including its path unrolled by five for unit strides. *)

  val dasum : int -> float array -> int -> int -> float
  (** [dasum n x xo incx] is the sum of the absolute values of the [n]
      elements of [x] from [xo] with stride [incx], and 0 when [n] is not
      positive or [incx] is not positive. The summation order is the
      Fortran's, including its path unrolled by six for a unit stride. *)

  val dnrm2 : int -> float array -> int -> int -> float
  (** [dnrm2 n x xo incx] is the Euclidean norm of the [n] element vector [x]
      from [xo] with stride [incx], computed by the scaled sum of squares so
      that it neither overflows nor underflows unnecessarily, and 0 when [n]
      is not positive. *)

  val idamax : int -> float array -> int -> int -> int
  (** [idamax n x xo incx] returns the index, counted from zero along the
      increment, of the first element of largest absolute value among the [n]
      elements of [x] starting at [xo] with stride [incx], and [-1] when [n]
      is below one or [incx] is not positive, which is where the Fortran
      returns its zero. *)
end

include module type of Raw
(** The Fortran shapes are in scope unqualified as well, because {!Lapack} and
    {!Linpack} and the differential drivers are all written that way. *)

(* ---- The OCaml surface ---- *)

val axpy : ?incx:int -> ?incy:int -> float -> float array -> float array
  -> unit
(** [axpy alpha x y] adds [alpha] times [x] to [y] in place. The increments
    default to one and have to be positive; at anything else they take every
    [inc] th element from the front, so the two vectors still need the same
    number of elements or you get [Invalid_argument]. *)

val dot : ?incx:int -> ?incy:int -> float array -> float array -> float
(** [dot x y] is the inner product. Increments read as in {!axpy}, and the
    summation order is the Fortran's. *)

val nrm2 : ?inc:int -> float array -> float
(** [nrm2 x] is the Euclidean norm, by the scaled sum of squares, so it
    doesn't overflow or underflow on its way to an answer that fits. [~inc]
    reads as in {!axpy}. *)

val asum : ?inc:int -> float array -> float
(** [asum x] is the sum of the absolute values, summed in the Fortran's
    order. [~inc] reads as in {!axpy}. *)

val scal : ?inc:int -> float -> float array -> unit
(** [scal alpha x] multiplies [x] by [alpha] in place. [~inc] reads as in
    {!axpy}, so with one it scales the lot and with two every other
    element. *)

val copy : ?incx:int -> ?incy:int -> float array -> float array -> unit
(** [copy x y] overwrites [y] with [x]. Increments read as in {!axpy}. If you
    just want a fresh array, [Array.copy] is the thing. *)

val swap : ?incx:int -> ?incy:int -> float array -> float array -> unit
(** [swap x y] exchanges the two in place. Increments read as in {!axpy}. *)

val rot : ?incx:int -> ?incy:int -> float array -> float array -> float
  -> float -> unit
(** [rot x y c s] applies the plane rotation with cosine [c] and sine [s] to
    the pairs drawn from [x] and [y], overwriting both. Increments read as in
    {!axpy}. *)

val rotg : float -> float -> float * float * float * float
(** [rotg a b] builds the plane rotation taking [(a, b)] onto the first axis
    and returns [(r, z, c, s)]: the rotated length, the single value [z] that
    [c] and [s] can be rebuilt from, then the cosine and the sine. *)

val iamax : ?inc:int -> float array -> int
(** [iamax x] is the zero based index, along [~inc], of the first element of
    largest absolute value. An empty vector has no such element, so that
    raises [Invalid_argument] rather than handing back the Fortran's
    sentinel. *)
