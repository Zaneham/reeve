(* Reeve, Copyright 2026 Zane Hambly.
   Polynomial interpolation, ported from the SLATEC double precision routines.
   A set of points goes through the divided difference table once and the
   table is then read to evaluate the polynomial, its derivatives, or its
   Taylor coefficients about a point. *)

module Raw : sig
  (** The Fortran shapes, argument for argument. These are what the
      differential in [diff/] compares against the reference Fortran. Use the
      functions below unless you specifically want to hold the table
      yourself. *)

  val dplint : int -> float array -> float array -> float array -> unit
  (** [dplint n x y c] fills the length [n] array [c] with the Newton divided
      difference table of the points [(x.(i), y.(i))], which {!dpolvl} and
      {!dpolcf} read afterwards. [x], [y] and [c] are length [n] and the
      abscissas must be distinct. Raises [Invalid_argument] for [n <= 0] and
      for a repeated abscissa, where the Fortran calls [XERMSG]. *)

  val dpolvl :
    int -> float -> float array -> int -> float array -> float array ->
    float array -> float * int
  (** [dpolvl nder xx yp n x c work] evaluates the polynomial that {!dplint}
      left in [c] at [xx], together with its first [nder] derivatives, which it
      writes to [yp] so that the derivative of order [j] lands in [yp.(j - 1)].
      [n], [x] and [c] must be as {!dplint} left them. [work] is scratch of
      length at least [2 * n] and is untouched when [nder <= 0]. Derivatives
      beyond the degree [n - 1] come back zero. Returns the pair
      [(yfit, ierr)], the value at [xx] and the Fortran error flag, which is 1
      for normal execution and is the only value the Fortran ever sets. *)

  val dpolcf :
    float -> int -> float array -> float array -> float array -> float array ->
    unit
  (** [dpolcf xx n x c d work] converts the table that {!dplint} left in [c]
      into the coefficients of the expansion about [xx], writing them to the
      length [n] array [d] so that
      [p z = d.(0) + d.(1) * (z - xx) + ... + d.(n - 1) * (z - xx) ** (n - 1)].
      [n], [x] and [c] must be as {!dplint} left them. [work] is scratch of
      length at least [2 * n]. {!dpolvl} is the more accurate way to evaluate
      the fit, since forming the Taylor coefficients loses accuracy. *)
end

(* ---- The OCaml surface ---- *)

type t
(** The polynomial through a set of points, held as its Newton divided
    difference table. You build one once and then evaluate it as often as you
    like. It keeps its own copy of the abscissas, so you can go on using the
    arrays you built it from. *)

val make : float array -> float array -> t
(** [make x y] is the polynomial of degree [Array.length x - 1] through the
    points [(x.(i), y.(i))]. The abscissas have to be distinct and there has
    to be at least one point, or you get [Invalid_argument]. *)

val points : t -> int
(** [points t] is how many points went into [t], so its degree is one less. *)

val eval : t -> float -> float
(** [eval t xx] is the value of [t] at [xx]. Nothing stops you evaluating well
    outside the points you fitted, and nothing good comes of it either. *)

val derivatives : t -> float -> int -> float * float array
(** [derivatives t xx nder] is the value of [t] at [xx] paired with its first
    [nder] derivatives, the one of order [j] at index [j - 1]. Derivatives past
    the degree come back zero, and a negative [nder] gets you
    [Invalid_argument]. *)

val coefficients : t -> float -> float array
(** [coefficients t xx] is the Taylor coefficients of [t] about [xx], so that
    [t] at [z] is [d.(0) + d.(1) * (z - xx) + ... + d.(n - 1) * (z - xx) ** (n - 1)].
    {!eval} is the more accurate way to get a value, because forming these
    loses accuracy. *)
