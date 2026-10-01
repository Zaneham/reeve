(* Reeve, Copyright 2026 Zane Hambly.
   Gamma, psi, the Airy functions and the Carlson elliptic integrals, ported
   from the SLATEC double precision routines. *)

val dai : float -> float
(** [dai x] is the Airy function Ai(x) for any finite [x]; it underflows to
    zero above x = 92.57 and raises nothing. *)

val dbi : float -> float
(** [dbi x] is the Bairy function Bi(x), the Airy function of the second kind.
    Raises [Invalid_argument] above x = 104.22, where Bi overflows. *)

val dcot : float -> float
(** [dcot x] is the cotangent of [x] in radians, by argument reduction onto
    the series for cot. Raises [Invalid_argument] where the result would
    overflow near zero, and where [x] is too large for the reduction to carry
    any precision. *)

val dcsevl : float -> float array -> int -> float
(** [dcsevl x cs n] evaluates the [n] term Chebyshev series [cs] at [x] by the
    backward recurrence, with [cs.(0)] the constant term. Raises
    [Invalid_argument] for [x] outside (-1,+1). *)

val dgamln : float -> float
(** [dgamln z] is the natural logarithm of the Gamma function for [z] > 0, by
    table lookup on the integers 1 to 100 and the asymptotic expansion
    elsewhere. Raises [Invalid_argument] for [z] <= 0. *)

val dpsi : float -> float
(** [dpsi x] is the digamma function psi(x). Raises [Invalid_argument] at x = 0
    and at the negative integers, and for large negative [x], where the
    cotangent term has no precision left. *)

val drc : float -> float -> float
(** [drc x y] is the degenerate Carlson elliptic integral RC(x,y), for [x] >= 0
    and [y] > 0. Raises [Invalid_argument] outside that domain, when
    max(x,y) exceeds max_float/5, when x+y falls below 5*min_float, or if the
    duplication iteration hits its bound. *)

val drd : float -> float -> float -> float
(** [drd x y z] is the Carlson elliptic integral of the second kind RD(x,y,z),
    for [x] >= 0, [y] >= 0, at most one of them zero, and [z] > 0. Raises
    [Invalid_argument] outside that domain, when max(x,y,z) exceeds roughly
    4.07e202, when min(x+y,z) falls below roughly 6.28e-206, or if the
    duplication iteration hits its bound. *)

val drf : float -> float -> float -> float
(** [drf x y z] is the Carlson elliptic integral of the first kind RF(x,y,z),
    for non-negative arguments of which at most one is zero. Raises
    [Invalid_argument] on a negative argument, when max(x,y,z) exceeds
    max_float/5, when min(x+y,x+z,y+z) falls below 5*min_float, or if the
    duplication iteration hits its bound. *)
