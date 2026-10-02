(* Reeve, Copyright 2026 Zane Hambly.
   Bessel functions,  that I've ported from the SLATEC double precision routines.
   Fixed order zero and one for I and K, and the Amos sequence routines
   for I, J and K of arbitrary non-negative order. *)

module Raw : sig
  (** The Fortran shapes of the sequence routines, argument for argument: you
      hand them the array and the integer [kode], and they hand you back the
      underflow count. These are what the differential in [diff/] compares
      against the reference Fortran. The rest of the module is already the
      shape an OCaml caller wants, so it isn't repeated here. *)

  val dbesi : float -> float -> int -> int -> float array -> int
  (** [dbesi x alpha kode n y] fills [y.(0) .. y.(n-1)] with the sequence
      I_(alpha+k)(x) for k = 0 .. n-1, unscaled for [kode = 1] and multiplied by
      [exp (-. x)] for [kode = 2], and returns the number of trailing components
      set to zero by underflow. Needs [x >= 0], [alpha >= 0], [n >= 1] and
      [kode] either 1 or 2; raises [Invalid_argument] otherwise, and also when
      [kode = 1] and [x] is large enough to overflow. *)

  val dbesj : float -> float -> int -> float array -> int
  (** [dbesj x alpha n y] fills [y.(0) .. y.(n-1)] with the sequence
      J_(alpha+k)(x) for k = 0 .. n-1, and returns the number of trailing
      components set to zero by underflow. Needs [x >= 0], [alpha >= 0] and
      [n >= 1]; raises [Invalid_argument] otherwise. *)

  val dbesk : float -> float -> int -> int -> float array -> int
  (** [dbesk x fnu kode n y] fills [y.(0) .. y.(n-1)] with the sequence
      K_(fnu+k)(x) for k = 0 .. n-1, unscaled for [kode = 1] and multiplied by
      [exp x] for [kode = 2], and returns the number of leading components set to
      zero by underflow. Needs [x > 0], [fnu >= 0], [n >= 1] and [kode] either 1
      or 2; raises [Invalid_argument] otherwise, and also when [fnu] or [n] is
      too large, or [x] too small, for the result to be representable. *)
end

(* ---- The OCaml surface ---- *)

val dbesi0 : float -> float
(** [dbesi0 x] is the modified Bessel function I_0(x) for any finite [x];
    raises [Invalid_argument] when [abs x] exceeds [log max_float]. *)

val dbesi0e : float -> float
(** [dbesi0e x] is the exponentially scaled I_0, [exp (-. abs x) *. i0 x], for
    any finite [x]; raises nothing. *)

val dbesi1 : float -> float
(** [dbesi1 x] is the modified Bessel function I_1(x) for any finite [x];
    raises [Invalid_argument] when [abs x] exceeds [log max_float]. *)

val dbesi1e : float -> float
(** [dbesi1e x] is the exponentially scaled I_1, [exp (-. abs x) *. i1 x], for
    any finite [x]; raises nothing. *)

val dbesk0 : float -> float
(** [dbesk0 x] is the modified Bessel function K_0(x) for [x > 0], returning
    zero once [x] is large enough for K_0 to underflow; raises
    [Invalid_argument] for [x <= 0]. *)

val dbesk0e : float -> float
(** [dbesk0e x] is the exponentially scaled K_0, [exp x *. k0 x], for [x > 0];
    raises [Invalid_argument] for [x <= 0]. *)

val dbesk1 : float -> float
(** [dbesk1 x] is the modified Bessel function K_1(x) for [x > 0], returning
    zero once [x] is large enough for K_1 to underflow; raises
    [Invalid_argument] for [x <= 0] or for [x] so small that K_1 overflows. *)

val dbesk1e : float -> float
(** [dbesk1e x] is the exponentially scaled K_1, [exp x *. k1 x], for [x > 0];
    raises [Invalid_argument] for [x <= 0] or for [x] so small that K_1
    overflows. *)

val dbesi : ?scaled:bool -> float -> float -> int -> float array * int
(** [dbesi x alpha n] is the [n] long sequence I_(alpha+k)(x) for
    k = 0 .. n-1, paired with the number of components at the end of it that
    underflowed to zero instead of being computed. [~scaled:true] multiplies
    the lot by [exp (-. x)], which is how you get an answer for an [x] that
    would otherwise overflow. Needs [x >= 0], [alpha >= 0] and [n >= 1];
    raises [Invalid_argument] otherwise, and on that overflow if you left the
    scaling off. *)

val dbesj : float -> float -> int -> float array * int
(** [dbesj x alpha n] is the [n] long sequence J_(alpha+k)(x) for
    k = 0 .. n-1, paired with the number of components at the end of it that
    underflowed to zero instead of being computed. Needs [x >= 0],
    [alpha >= 0] and [n >= 1]; raises [Invalid_argument] otherwise. *)

val dbesk : ?scaled:bool -> float -> float -> int -> float array * int
(** [dbesk x fnu n] is the [n] long sequence K_(fnu+k)(x) for k = 0 .. n-1,
    paired with the number of components that underflowed to zero. They're the
    ones at the start this time, since K grows with its order and it's the low
    orders that go under. [~scaled:true] multiplies the lot by [exp x]. Needs
    [x > 0], [fnu >= 0] and [n >= 1]; raises [Invalid_argument] otherwise, and
    when [fnu] or [n] is too large, or [x] too small, for the result to be
    representable. *)
