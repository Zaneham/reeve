(* Reeve, Copyright 2026 Zane Hambly.
   Quadrature, ported from the SLATEC double precision routines. Where the
   Fortran took an EXTERNAL integrand these take a [float -> float] closure.
   A matrix of Taylor derivatives is a [float array] in column major order
   with an explicit leading dimension, so the element (i, j) zero based lives
   at [i + j * ldc]. *)

module Raw : sig
  (** The Fortran shapes, argument for argument. Every routine here reports its
      status through the returned [ierr], keeping the Fortran's codes, and
      raises only where the Fortran calls the input improper. The two
      piecewise polynomial routines are only here, because a piecewise
      polynomial is still a bare array with a leading dimension. *)

  val davint : float array -> float array -> int -> float -> float
    -> float * int
  (** [davint x y n xlo xup] integrates the function tabulated as the [n]
      abscissas [x], which must increase, against the values [y], over [xlo] to
      [xup], by overlapping parabolas, or by the trapezoid rule when [n] is 2.
      Returns [(ans, ierr)]. [ierr] is 1 when the integration was performed, 2
      when [xup] is below [xlo], 3 when fewer than three abscissas lie between
      the limits, 4 when [x] does not strictly increase and 5 when [n] is below
      2, and [ans] is zero for every code but 1. *)

  val dgaus8 : (float -> float) -> float -> float -> float
    -> float * int * float
  (** [dgaus8 f a b err] integrates [f] from [a] to [b], which may be below
      [a], by an adaptive eight point Legendre-Gauss rule to the pseudorelative
      tolerance [err]. Returns [(ans, ierr, err_out)]. [ierr] is 1 when [ans]
      most likely meets the tolerance or [a] equals [b], -1 when [a] and [b]
      are too nearly equal to integrate and [ans] is zero, and 2 when [ans]
      probably does not meet the tolerance. [err_out] is an estimate of the
      absolute error in [ans] when [err] was negative, and [err] itself
      otherwise. *)

  val dqnc79 : (float -> float) -> float -> float -> float
    -> float * int * int
  (** [dqnc79 f a b err] integrates [f] from [a] to [b], which may be below
      [a], by an adaptive seven point Newton-Cotes rule to the tolerance [err].
      Returns [(ans, ierr, k)], where [k] counts the evaluations of [f] and
      [ierr] carries the codes of {!dgaus8}, with -1 also covering [a] equal to
      [b]. *)

  val dppgq8 :
    (float -> float) -> int -> float array -> float array -> int -> int -> int
    -> float -> float -> float -> float * int * float
  (** [dppgq8 f ldc c xi lxi kk id a b err] integrates the product of [f] and
      the [id]th derivative of the piecewise polynomial of order [kk] held as
      the right Taylor derivatives [c], of shape [kk] by [lxi] with leading
      dimension [ldc], over the [lxi] pieces broken at the [lxi + 1] points
      [xi], from [a] to [b]. Returns [(ans, ierr, err_out)] with the meanings
      of {!dgaus8}. Raises [Invalid_argument] unless [kk] is at least 1, [ldc]
      at least [kk], [lxi] at least 1 and [0 <= id < kk]. *)

  val dpfqad :
    (float -> float) -> int -> float array -> float array -> int -> int -> int
    -> float -> float -> float -> float * int
  (** [dpfqad f ldc c xi lxi k id x1 x2 tol] integrates the product of [f] and
      the [id]th derivative of the piecewise polynomial [(c, xi, lxi, k)],
      described as in {!dppgq8}, from [x1] to [x2], splitting the interval at
      the breakpoints it contains and handing each piece to {!dppgq8}. Returns
      [(quad, ierr)], where [ierr] is 1 normally and 2 when some piece did not
      meet [tol]. Raises [Invalid_argument] unless [k] is at least 1, [ldc] at
      least [k], [lxi] at least 1, [0 <= id < k], and [tol] lies between the
      unit roundoff and 0.1. *)
end

(* ---- The OCaml surface ---- *)

type result = { value : float; error : float option; converged : bool }
(** What an adaptive rule has to say. [value] is the integral. [converged] is
    false when the rule ran out of subdivisions before it met the tolerance, in
    which case [value] is still the best it has, just not as good as you asked
    for. [error] is an estimate of the absolute error, and is [None] for a rule
    that doesn't make one. *)

exception Too_narrow of float * float
(** The two limits are so nearly equal that the rule can't tell them apart
    well enough to integrate between them. They aren't equal, though: an
    interval of no width integrates to zero. *)

val gauss8 : ?tol:float -> (float -> float) -> float -> float -> result
(** [gauss8 f a b] integrates [f] from [a] to [b], either way round, by an
    adaptive eight point Legendre-Gauss rule. [~tol] is pseudorelative and
    defaults to about the square root of the unit roundoff. Raises
    {!Too_narrow}. *)

val newton_cotes7 : ?tol:float -> (float -> float) -> float -> float -> result
(** [newton_cotes7 f a b] is {!gauss8}'s job done by an adaptive seven point
    Newton-Cotes rule instead, which tends to want more evaluations of [f] but
    copes better with an integrand that isn't smooth. It gives no error
    estimate. Raises {!Too_narrow}. *)

val integrate_table : float array -> float array -> float -> float -> float
(** [integrate_table x y lo up] integrates the function tabulated as the points
    [(x.(i), y.(i))] over [lo] to [up] by overlapping parabolas, or by the
    trapezoid rule if you only hand it two points. The abscissas have to
    increase, and there have to be three of them between the limits unless
    there are only two points in all, or you get [Invalid_argument]. *)
