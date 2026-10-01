(* Reeve, Copyright 2026 Zane Hambly.
   The exponential integrals E1 and Ei, the logarithmic integral, Spence's
   integral, Dawson's function, and the elementary kernels SLATEC carries
   alongside them, ported from the SLATEC double precision routines. *)

val de1 : float -> float
(** [de1 x] is the exponential integral E1(x) for [x] > 0 and its Cauchy
    principal value for [x] < 0, returning zero once [x] is large enough for
    E1 to underflow. Raises [Invalid_argument] at x = 0. *)

val dei : float -> float
(** [dei x] is the exponential integral Ei(x) for [x] > 0 and its Cauchy
    principal value for [x] < 0, computed as [-. de1 (-. x)]. Raises
    [Invalid_argument] at x = 0. *)

val dexint : float -> int -> int -> int -> float -> float array -> int
(** [dexint x n kode m tol en] fills [en.(0) .. en.(m-1)] with the sequence
    E(n+k, x) for k = 0 .. m-1, unscaled for [kode = 1] and multiplied by
    [exp x] for [kode = 2], and returns the number of components set to zero
    by underflow, which is either zero or [m]. Needs [x >= 0], [n >= 1],
    [m >= 1], [kode] either 1 or 2, [tol] between the double precision unit
    roundoff and 0.1, and not both [x] = 0 and [n] = 1; raises
    [Invalid_argument] otherwise, and also where the Fortran reports its
    algorithm termination condition unmet. *)

val dli : float -> float
(** [dli x] is the logarithmic integral li(x) for [x] > 0, computed as
    [dei (log x)]. Raises [Invalid_argument] for [x] <= 0 and at x = 1, where
    li has its pole. *)

val dspenc : float -> float
(** [dspenc x] is Spence's integral in K. Mitchell's form, the integral of
    -log(1-y)/y from 0 to [x], for any finite [x]; raises nothing. *)

val ddaws : float -> float
(** [ddaws x] is Dawson's integral, exp(-x^2) times the integral of exp(t^2)
    from 0 to [x], for any finite [x], returning zero once [abs x] is large
    enough for the result to underflow; raises nothing. *)

val dpsixn : int -> float
(** [dpsixn n] is psi(n) for the integer [n] >= 1, by table lookup up to 100
    and from the asymptotic expansion above that. It is SLATEC's subsidiary to
    {!dexint}. Raises [Invalid_argument] for [n] < 1. *)

val dexprl : float -> float
(** [dexprl x] is the relative error exponential (exp x - 1)/x, by Taylor
    series near zero and directly elsewhere; raises nothing. *)

val dlnrel : float -> float
(** [dlnrel x] is log(1+x), accurate in the sense of relative error near
    x = 0. Raises [Invalid_argument] for [x] <= -1. *)

val dcbrt : float -> float
(** [dcbrt x] is the cube root of [x], by a single precision polynomial
    starting value refined by Newton iteration; raises nothing. *)

val dsindg : float -> float
(** [dsindg x] is the sine of [x] measured in degrees, exact at the multiples
    of 90; raises nothing. *)

val dcosdg : float -> float
(** [dcosdg x] is the cosine of [x] measured in degrees, exact at the
    multiples of 90; raises nothing. *)

val d9atn1 : float -> float
(** [d9atn1 x] is (atan x - x)/x^3, the arc tangent from first order with
    relative accuracy, so that atan x = x + x^3 * d9atn1 x. Raises
    [Invalid_argument] when [abs x] exceeds 1.571/eps, where the answer has no
    precision left. *)

val d9ln2r : float -> float
(** [d9ln2r x] is (log(1+x) - x + x^2/2)/x^3, the logarithm from second order
    with relative accuracy, so that log(1+x) = x - x^2/2 + x^3 * d9ln2r x.
    Raises [Invalid_argument] for [x] above roughly 7.2e8, where the answer has
    no precision left. *)

val d9pak : float -> int -> float
(** [d9pak y n] is [y] scaled by two to the power [n], returning zero where
    that would underflow below the smallest normal. Raises [Invalid_argument]
    where it would overflow. *)

val d9upak : float -> float * int
(** [d9upak x] is the pair [(y, n)] with [x = y * 2 ** n] and [abs y] in
    [[0.5, 1.0)], or [(x, 0)] for a zero, an infinity or a NaN. *)
