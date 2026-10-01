(* Reeve, Copyright 2026 Zane Hambly.
   Extended-range arithmetic of Smith, Olver and Lozier, and the Legendre
   functions of the first and second kind built on top of it, ported from the
   SLATEC double precision routines. *)

type xnum = { x : float; ix : int }
(** An extended-range number [{ x; ix }] stands for the real number
    [x *. radix ** ix], where [radix] is the radix of the arithmetic the
    enclosing {!env} was set up for. Every real number has many such forms;
    [dxadj] puts one in adjusted form and [dxred] recovers an ordinary float
    where that is possible. Only [dxcon] departs from this reading: its result
    is to be read in base ten. *)

type env
(** The machine-dependent constants the extended-range routines need: radix,
    overflow limits, the index bound and the base ten conversion table. It can
    only be made by {!dxset}, so no routine here can be called before the
    package is initialised. *)

type kind = Pneg | Q | Ppos | Pnorm
(** Which Legendre function {!dxlegf} is to compute, the Fortran [ID]:
    [Pneg] is P of negative order P(-mu,nu,x) (ID 1), [Q] is the function of
    the second kind Q(mu,nu,x) (ID 2), [Ppos] is P of positive order
    P(mu,nu,x) (ID 3) and [Pnorm] is the normalised Legendre polynomial
    (ID 4). *)

type mode = By_x | By_theta
(** How {!dxnrmp} reads its argument, the Fortran [MODE]: [By_x] takes it as
    x itself (MODE 1), [By_theta] as an angle whose cosine is x (MODE 2). *)

val dxset : int -> int -> float -> int -> env
(** [dxset irad nradpl dzero nbits] builds the extended-range parameters.
    [irad] is the radix of the floating-point arithmetic, [nradpl] the number
    of radix places it carries, [dzero] the smallest of 1/min_float,
    max_float and the largest float whose log10 can be formed, and [nbits] the
    number of bits excluding the sign in an integer. Any argument given as 0
    (0.0 for [dzero]) is taken from the machine instead, which is the normal
    call: [dxset 0 0 0.0 0]. Raises [Invalid_argument] for a radix that is not
    2, 4, 8 or 16, for [nbits] outside 15 to {!Sys.int_size} - 1, for
    [nradpl] outside 1 to 120 radix-2 places or not less than the resulting
    scaling exponent, for a [dzero] that leaves too little exponent range, and
    where radix raised to twice that exponent would overflow. *)

val dxadj : env -> xnum -> xnum
(** [dxadj env a] returns [a] in adjusted form, so that its principal part is
    zero or satisfies radix**(-l) <= |x| < radix**l, where l is the scaling
    exponent [dxset] chose. Two adjusted numbers can be multiplied or divided
    without the principal part overflowing. The principal part of [a] must
    already lie between radix**(-2l) and radix**(2l). Raises
    [Invalid_argument] if the auxiliary index would leave its range. *)

val dxadd : env -> xnum -> xnum -> xnum
(** [dxadd env a b] is the extended-range sum of [a] and [b], returned in
    adjusted form. The operands need not be adjusted, but their principal
    parts must lie between radix**(-2l) and radix**(2l). Raises
    [Invalid_argument] if the auxiliary index of the sum would leave its
    range. *)

val dxred : env -> xnum -> xnum
(** [dxred env a] returns [a] with its auxiliary index reduced to zero, so
    that the whole value sits in the principal part and ordinary float
    arithmetic can carry on with it. If [a] lies outside radix**(-2l) to
    radix**(2l) it is returned untouched, index and all. Raises
    [Invalid_argument] if a scaling loop hits its bound, which a principal
    part that is not finite will do when the index is nonzero. *)

val dxcon : env -> xnum -> xnum
(** [dxcon env a] converts [a] to the decimal extended-range form, ready for
    printing: the result is to be read as [x *. 10.0 ** ix] with
    1/10 <= |x| < 1, except that a value within radix**(-2l) to radix**(2l)
    comes back reduced, with index zero. Raises [Invalid_argument] if the
    index is too large for the conversion table or a scaling loop hits its
    bound. *)

val dxlegf : env -> float -> int -> int -> int -> float -> kind -> xnum array
(** [dxlegf env dnu1 nudiff mu1 mu2 theta kind] computes a vector of Legendre
    functions of [kind] at x = cos [theta]. With [nudiff] = 0 it runs over the
    orders mu = [mu1], [mu1]+1, ..., [mu2] at degree nu = [dnu1]; with
    [mu1] = [mu2] it runs over the degrees nu = [dnu1], [dnu1]+1, ...,
    [dnu1]+[nudiff] at order mu = [mu1]. The result has length
    [mu2]-[mu1]+[nudiff]+1 and is in reduced form wherever the value fits a
    float, so an element whose index is zero needs no further thought about
    extended range. Requires [dnu1] >= -0.5, [nudiff] >= 0,
    0 <= [mu1] <= [mu2], [theta] in (0, pi/2], and either [nudiff] = 0 or
    [mu1] = [mu2]; [Pnorm] additionally requires an integer [dnu1]. Raises
    [Invalid_argument] otherwise, and for an extended-range index overflow
    underneath. *)

val dxnrmp : env -> int -> int -> int -> float -> mode -> xnum array * int
(** [dxnrmp env nu mu1 mu2 darg mode] computes the normalised Legendre
    polynomials of degree [nu] and orders [mu1] to [mu2] at the argument
    [darg], read according to [mode]. The normalisation is the one that makes
    the integral of the square over -1 to 1 equal to one. Returns the vector,
    of length [mu2]-[mu1]+1 and in reduced form wherever the value fits a
    float, together with an estimate of the number of decimal digits lost to
    rounding, so a result good to d digits in [darg] is good to d minus that
    many. Requires [nu] >= 0, 0 <= [mu1] <= [mu2], and |[darg]| <= 1 for
    [By_x] or |[darg]| <= pi for [By_theta]. Raises [Invalid_argument]
    otherwise. *)
