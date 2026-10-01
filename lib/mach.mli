(* Reeve, Copyright 2026 Zane Hambly.
   The machine constants the ported routines ask D1MACH and DLAMCH for, the
   float32 rounding the Fortran gets from a single precision expression, and
   the Blue scaling thresholds that follow from the exponent range. This file
   depends on nothing. *)

val eps_dp : float
(** The relative spacing one above, D1MACH(4), which is [epsilon_float]. *)

val eps_2_dp : float
(** Half of {!eps_dp}, D1MACH(3), the largest relative error of a correctly
    rounded result. The SLATEC series truncation thresholds are stated in
    terms of this one. *)

val tiny_dp : float
(** The smallest positive normal double, D1MACH(1). *)

val huge_dp : float
(** The largest finite double, D1MACH(2). *)

val log10_radix_dp : float
(** The base ten logarithm of the radix, D1MACH(5), so [log10 2.0]. *)

val digits_dp : int
(** The number of base two digits in the significand, I1MACH(14), so 53. *)

val min_exp_dp : int
(** The smallest exponent of the radix giving a normal number, I1MACH(15),
    so -1021. *)

val max_exp_dp : int
(** The largest exponent of the radix giving a finite number, I1MACH(16),
    so 1024. *)

val single : float -> float
(** [single x] is [x] rounded to the nearest float32 and widened back, which
    is what a single precision subexpression in the Fortran computes. The
    SLATEC term counts and a few thresholds are written in single precision
    and the rounding is visible in the result. *)

val fdigits : int
(** The significand digit count as DLAMCH('N') derives it, from [frexp] of
    the epsilon rather than from a literal. *)

val fminexp : int
(** The minimum exponent as DLAMCH('M') derives it, from [frexp] of the
    smallest normal. *)

val fmaxexp : int
(** The maximum exponent as DLAMCH('L') derives it, from [frexp] of the
    largest finite double. *)

val tsml : float
(** The lower threshold of Blue's three range sum of squares: below this an
    element is accumulated scaled up by {!ssml}. *)

val tbig : float
(** The upper threshold of Blue's three range sum of squares: above this an
    element is accumulated scaled down by {!sbig}. *)

val ssml : float
(** The scaling applied to an element below {!tsml}, a power of two, so the
    scaling is exact. *)

val sbig : float
(** The scaling applied to an element above {!tbig}, a power of two, so the
    scaling is exact. *)
