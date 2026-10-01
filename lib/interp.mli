(* Reeve, Copyright 2026 Zane Hambly.
   Polynomial interpolation, ported from the SLATEC double precision routines.
   A set of points goes through {!dplint} once and the table it leaves in [c]
   is then read by {!dpolvl} and {!dpolcf}. *)

val dplint : int -> float array -> float array -> float array -> unit
(** [dplint n x y c] fills the length [n] array [c] with the Newton divided
    difference table of the points [(x.(i), y.(i))], which {!dpolvl} and
    {!dpolcf} read afterwards. [x], [y] and [c] are length [n] and the
    abscissas must be distinct. Raises [Invalid_argument] for [n <= 0] and for
    a repeated abscissa, where the Fortran calls [XERMSG]. *)

val dpolvl :
  int -> float -> float array -> int -> float array -> float array ->
  float array -> float * int
(** [dpolvl nder xx yp n x c work] evaluates the polynomial that {!dplint} left
    in [c] at [xx], together with its first [nder] derivatives, which it writes
    to [yp] so that the derivative of order [j] lands in [yp.(j - 1)]. [n], [x]
    and [c] must be as {!dplint} left them. [work] is scratch of length at
    least [2 * n] and is untouched when [nder <= 0]. Derivatives beyond the
    degree [n - 1] come back zero. Returns the pair [(yfit, ierr)], the value
    at [xx] and the Fortran error flag, which is 1 for normal execution and is
    the only value the Fortran ever sets. *)

val dpolcf :
  float -> int -> float array -> float array -> float array -> float array ->
  unit
(** [dpolcf xx n x c d work] converts the table that {!dplint} left in [c] into
    the coefficients of the expansion about [xx], writing them to the length
    [n] array [d] so that
    [p z = d.(0) + d.(1) * (z - xx) + ... + d.(n - 1) * (z - xx) ** (n - 1)].
    [n], [x] and [c] must be as {!dplint} left them. [work] is scratch of
    length at least [2 * n]. {!dpolvl} is the more accurate way to evaluate the
    fit, since forming the Taylor coefficients loses accuracy. *)
