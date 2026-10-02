(* Reeve, Copyright 2026 Zane Hambly.
   The gamma and beta family: the reciprocal and log gamma, factorials and
   binomials, Pochhammer's symbol, the complete and incomplete beta, the
   incomplete gamma trio, the polygamma sequence, Kummer's confluent
   hypergeometric function U, and the Wigner 3j and 6j coefficients, ported
   from the SLATEC double precision routines. *)

module Raw : sig
  (** The Fortran shapes of the four routines that fill an array you pass them,
      argument for argument. These are what the differential in [diff/]
      compares against the reference Fortran. The rest of the module is already
      the shape an OCaml caller wants, so it isn't repeated here. *)

  val dpsifn : float -> int -> int -> int -> float array -> int
  (** [dpsifn x n kode m ans] fills [ans.(0) .. ans.(m-1)] with the scaled
      derivatives ((-1)^(k+1)/Gamma(k+1)) psi(k,x) for k = n .. n+m-1, except
      that for [kode = 2] and [n = 0] the first component is -psi(x) + log x
      instead of -psi(x), and returns the number of trailing components set to
      zero by underflow. Needs [x] > 0, [n] >= 0, [m] >= 1, [kode] either 1 or
      2 and [ans] at least [m] long; raises [Invalid_argument] otherwise, and
      when [x] is too small or [n+m-1] too large for the result to be
      representable. *)

  val drc3jj :
    float -> float -> float -> float -> float array -> int -> float * float
  (** [drc3jj l2 l3 m2 m3 thrcof ndim] fills [thrcof] with the Wigner 3j symbols
      (l1 l2 l3 / -m2-m3 m2 m3) for l1 running over its allowed values in unit
      steps, and returns the pair (l1min, l1max) of those limits, so that
      [thrcof.(0)] is the symbol at l1min and l1max - l1min + 1 components are
      set. Needs [ndim] at least that many, [ndim] no more than the length of
      [thrcof], l2 >= |m2| and l3 >= |m3| with l2+|m2| and l3+|m3| integral;
      raises [Invalid_argument] otherwise. *)

  val drc3jm :
    float -> float -> float -> float -> float array -> int -> float * float
  (** [drc3jm l1 l2 l3 m1 thrcof ndim] fills [thrcof] with the Wigner 3j symbols
      (l1 l2 l3 / m1 m2 -m1-m2) for m2 running over its allowed values in unit
      steps, and returns the pair (m2min, m2max) of those limits, so that
      [thrcof.(0)] is the symbol at m2min and m2max - m2min + 1 components are
      set. Needs [ndim] at least that many, [ndim] no more than the length of
      [thrcof], l1 >= |m1| with l1+|m1| integral, l1, l2 and l3 triangular and
      l1+l2+l3 integral; raises [Invalid_argument] otherwise. *)

  val drc6j :
    float ->
    float ->
    float ->
    float ->
    float ->
    float array ->
    int ->
    float * float
  (** [drc6j l2 l3 l4 l5 l6 sixcof ndim] fills [sixcof] with the Wigner 6j
      symbols (l1 l2 l3 / l4 l5 l6) for l1 running over its allowed values in
      unit steps, and returns the pair (l1min, l1max) of those limits, so that
      [sixcof.(0)] is the symbol at l1min and l1max - l1min + 1 components are
      set. Needs [ndim] at least that many, [ndim] no more than the length of
      [sixcof], the triads (l4,l2,l6) and (l4,l5,l3) triangular and l2+l3+l5+l6
      and l4+l2+l6 integral; raises [Invalid_argument] otherwise. *)
end

(* ---- The OCaml surface ---- *)

val dgamlm : unit -> float * float
(** [dgamlm ()] is the pair (xmin, xmax) bounding the argument of the Gamma
    function, xmin the smallest value below which Gamma underflows and xmax the
    largest above which it overflows. Raises [Invalid_argument] if either
    Newton iteration fails to settle in ten steps. *)

val dgamr : float -> float
(** [dgamr x] is the reciprocal 1/Gamma(x), returning zero at zero and at the
    negative integers where Gamma has its poles. Raises nothing. *)

val dlgams : float -> float * float
(** [dlgams x] is the pair (log |Gamma(x)|, sign of Gamma(x)), the sign being
    +1.0 or -1.0. Raises [Invalid_argument] at zero and at the negative
    integers. *)

val d9lgmc : float -> float
(** [d9lgmc x] is the log gamma correction term for [x] >= 10, that is
    log Gamma(x) - (x - 0.5) log x + x - 0.5 log (2 pi). Raises
    [Invalid_argument] for [x] < 10. *)

val dfac : int -> float
(** [dfac n] is the factorial of [n], exact from a table up to 30 and from the
    Stirling series above that. Raises [Invalid_argument] for [n] < 0 and for
    [n] > 170, where the factorial overflows. *)

val dbinom : int -> int -> float
(** [dbinom n m] is the binomial coefficient n choose m, rounded to an integer
    whenever the result is small enough to be one exactly. Needs
    [0 <= m <= n]; raises [Invalid_argument] otherwise, and when the result
    overflows. *)

val dpoch : float -> float -> float
(** [dpoch a x] is Pochhammer's generalised symbol Gamma(a+x)/Gamma(a),
    returning zero when [a] is a non-positive integer and [a+x] is not. Raises
    [Invalid_argument] when [a+x] is a non-positive integer and [a] is not, so
    that the ratio has a pole. *)

val dpoch1 : float -> float -> float
(** [dpoch1 a x] is (Gamma(a+x)/Gamma(a) - 1)/x, which at [x = 0] is psi(a),
    computed so as to keep the accuracy that the plain difference loses for
    small [x]. Raises [Invalid_argument] where [dpoch] or [dpsi] does, and if
    the Bernoulli series would need more than twenty terms. *)

val dlbeta : float -> float -> float
(** [dlbeta a b] is the natural logarithm of the complete beta function
    B(a,b), for [a] > 0 and [b] > 0. Raises [Invalid_argument] otherwise. *)

val dbeta : float -> float -> float
(** [dbeta a b] is the complete beta function B(a,b), for [a] > 0 and [b] > 0,
    underflowing to zero once the logarithm falls below log min_float. Raises
    [Invalid_argument] for a non-positive argument. *)

val dbetai : float -> float -> float -> float
(** [dbetai x p q] is the incomplete beta function ratio I_x(p,q), the integral
    of t^(p-1) (1-t)^(q-1) from zero to [x] divided by B(p,q), clamped to the
    interval from zero to one. Needs [0 <= x <= 1], [p] > 0 and [q] > 0; raises
    [Invalid_argument] otherwise. *)

val dgami : float -> float -> float
(** [dgami a x] is the incomplete gamma function, the integral of
    t^(a-1) exp(-t) from zero to [x], for [a] > 0 and [x] >= 0. Raises
    [Invalid_argument] outside that domain. *)

val dgamic : float -> float -> float
(** [dgamic a x] is the complementary incomplete gamma function, the integral
    of t^(a-1) exp(-t) from [x] to infinity, for [x] >= 0 and any [a] when
    [x] > 0. Raises [Invalid_argument] for negative [x], and for [x = 0] with
    [a] <= 0 where the integral diverges. *)

val dgamit : float -> float -> float
(** [dgamit a x] is Tricomi's form of the incomplete gamma function,
    x^(-a) P(a,x) continued analytically in [a], which stays finite where
    [dgami] does not. Needs [x] >= 0; raises [Invalid_argument] otherwise. *)

val d9gmic : float -> float -> float -> float
(** [d9gmic a x alx] is the complementary incomplete gamma function for [a]
    near a negative integer and [x] small, with [alx] the logarithm of [x].
    Raises [Invalid_argument] for [a] > 0, for [x] <= 0, and if the series
    fails to converge in 200 terms. *)

val d9gmit : float -> float -> float -> float -> float
(** [d9gmit a x algap1 sgngam] is Tricomi's incomplete gamma function for
    [x] <= 1, taking [algap1] as log |Gamma(a+1)| and [sgngam] as the sign of
    Gamma(a+1). Raises [Invalid_argument] for [x] <= 0 and if the Taylor series
    fails to converge in 200 terms. *)

val d9lgic : float -> float -> float -> float
(** [d9lgic a x alx] is the logarithm of the complementary incomplete gamma
    function for large [x] and for [a] <= [x], with [alx] the logarithm of
    [x]. Raises [Invalid_argument] if the continued fraction fails to converge
    in 300 terms. *)

val d9lgit : float -> float -> float -> float
(** [d9lgit a x algap1] is the logarithm of Tricomi's incomplete gamma function
    for perfect precision, taking [algap1] as log Gamma(a+1). Needs [x] > 0 and
    [x] <= [a]; raises [Invalid_argument] otherwise and if the continued
    fraction fails to converge in 200 terms. *)

val dpsixn : int -> float
(** [dpsixn n] is psi(n) for the integer [n] >= 1, by table lookup up to 100
    and from the asymptotic expansion above that. Raises [Invalid_argument] for
    [n] < 1. *)

val dpsifn : ?scaled:bool -> float -> int -> int -> float array * int
(** [dpsifn x n m] is the [m] long sequence of scaled derivatives
    ((-1)^(k+1)/Gamma(k+1)) psi(k,x) for k = n .. n+m-1, paired with the number
    of components at the end of it that underflowed to zero instead of being
    computed. With [~scaled:true] and [n = 0] the first component comes back as
    -psi(x) + log x rather than -psi(x). Needs [x] > 0, [n] >= 0 and [m] >= 1;
    raises [Invalid_argument] otherwise, and when [x] is too small or [n+m-1]
    too large for the result to be representable. *)

val dchu : float -> float -> float -> float
(** [dchu a b x] is Kummer's confluent hypergeometric function U(a,b,x) of the
    second kind, for [x] > 0. Raises [Invalid_argument] for [x] <= 0, when
    1+a-b is so near zero that the algorithm loses its footing for small [x],
    and if either ascending series fails to converge in 1000 terms. *)

val d9chu : float -> float -> float -> float
(** [d9chu a b z] is the logarithmic-derivative-free continued fraction for
    U(a,b,z) scaled by z^a, that is z^a U(a,b,z), used by [dchu] for large [z].
    Raises [Invalid_argument] if the fraction fails to converge in 300
    terms. *)

val drc3jj : float -> float -> float -> float -> float array * float * float
(** [drc3jj l2 l3 m2 m3] is the triple (coefficients, l1min, l1max): the Wigner
    3j symbols (l1 l2 l3 / -m2-m3 m2 m3) for l1 running from l1min to l1max in
    unit steps, so the array is l1max - l1min + 1 long and its first element is
    the symbol at l1min. Needs l2 >= |m2| and l3 >= |m3| with l2+|m2| and
    l3+|m3| integral; raises [Invalid_argument] otherwise, and if the run of
    coefficients is too long to allocate. *)

val drc3jm : float -> float -> float -> float -> float array * float * float
(** [drc3jm l1 l2 l3 m1] is the triple (coefficients, m2min, m2max): the Wigner
    3j symbols (l1 l2 l3 / m1 m2 -m1-m2) for m2 running from m2min to m2max in
    unit steps, so the array is m2max - m2min + 1 long and its first element is
    the symbol at m2min. Needs l1 >= |m1| with l1+|m1| integral, l1, l2 and l3
    triangular and l1+l2+l3 integral; raises [Invalid_argument] otherwise, and
    if the run of coefficients is too long to allocate. *)

val drc6j :
  float -> float -> float -> float -> float -> float array * float * float
(** [drc6j l2 l3 l4 l5 l6] is the triple (coefficients, l1min, l1max): the
    Wigner 6j symbols (l1 l2 l3 / l4 l5 l6) for l1 running from l1min to l1max
    in unit steps, so the array is l1max - l1min + 1 long and its first element
    is the symbol at l1min. Needs the triads (l4,l2,l6) and (l4,l5,l3)
    triangular and l2+l3+l5+l6 and l4+l2+l6 integral; raises
    [Invalid_argument] otherwise, and if the run of coefficients is too long to
    allocate. *)
