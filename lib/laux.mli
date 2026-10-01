(* Reeve, Copyright 2026 Zane Hambly.
   The LAPACK auxiliary routines, ported from the reference double precision
   Fortran. A matrix is a [float array] in column major order with a base
   offset [ao] and an explicit leading dimension [lda], so the element (i, j)
   zero based lives at [ao + i + j * lda]. Row and pivot indices are zero
   based throughout. *)

type machine_parameter =
  | Eps
  | Safe_min
  | Base
  | Precision
  | Digits
  | Round
  | Min_exponent
  | Underflow
  | Max_exponent
  | Overflow
(** The machine parameter selected by {!dlamch}, in place of the Fortran
    character: [Eps] for 'E', [Safe_min] for 'S', [Base] for 'B',
    [Precision] for 'P', [Digits] for 'N', [Round] for 'R',
    [Min_exponent] for 'M', [Underflow] for 'U', [Max_exponent] for 'L' and
    [Overflow] for 'O'. *)

type uplo = Upper | Lower | Full
(** Which part of a rectangular array a routine touches: [Upper] and [Lower]
    for the triangular or trapezoidal part, [Full] for the whole leading
    m by n submatrix, in place of the Fortran [uplo] character. *)

type norm = Max_abs | One_norm | Inf_norm | Frobenius
(** The norm selected by {!dlange}, in place of the Fortran character:
    [Max_abs] for 'M', [One_norm] for '1' or 'O', [Inf_norm] for 'I' and
    [Frobenius] for 'F' or 'E'. *)

val xerbla : string -> int -> 'a
(** [xerbla srname info] is the LAPACK error reporter. Where the Fortran
    prints and stops, this raises [Invalid_argument] naming the routine
    [srname] and the position [info] of the offending argument, counting from
    one as the Fortran does. *)

val dlamch : machine_parameter -> float
(** [dlamch cmach] returns the selected IEEE double parameter: the relative
    machine precision, the safe minimum whose reciprocal does not overflow,
    the base of the arithmetic, [Eps * Base], the number of significant base
    digits, 1 for round rather than chop, the minimum and maximum exponents
    before gradual underflow and before overflow, and the underflow and
    overflow thresholds. *)

val dlaisnan : float -> float -> bool
(** [dlaisnan din1 din2] is the inequality [din1 <> din2], which is true for
    any NaN operand and so detects a NaN when the same value is passed
    twice. *)

val disnan : float -> bool
(** [disnan din] is true when [din] is NaN. *)

val dlapy2 : float -> float -> float
(** [dlapy2 x y] returns [sqrt (x * x + y * y)] without unnecessary overflow
    or underflow, and returns the NaN operand when either is NaN. *)

val dlapy3 : float -> float -> float -> float
(** [dlapy3 x y z] returns [sqrt (x * x + y * y + z * z)] without unnecessary
    overflow or underflow. *)

val dlassq : int -> float array -> int -> int -> float -> float -> float * float
(** [dlassq n x xo incx scale sumsq] updates the scaled sum of squares of the
    [n] element vector in [x] from base offset [xo] with increment [incx],
    and returns the new [(scale, sumsq)] such that
    [scale * scale * sumsq] is the sum of squares of the vector plus
    [scale * scale * sumsq] on entry. Returns the arguments unchanged when
    either [scale] or [sumsq] is NaN. *)

val dlaswp : int -> float array -> int -> int -> int -> int -> int array
  -> int -> int -> unit
(** [dlaswp n a ao lda k1 k2 ipiv ipo incx] applies a series of row
    interchanges to the [n] columns of [a], exchanging row [i] with row
    [ipiv (ipo + k1 + (i - k1) * abs incx)] for each row [i] from [k1] to
    [k2] inclusive, all indices zero based. [ipo] is the base offset of the
    pivot vector, which a recursive factorisation shifts as it descends.
    [incx] of 1 walks [ipiv] forward from [k1] and -1 walks it backward,
    applying the interchanges in reverse order, while 0 does nothing. *)

val dlaset : uplo -> int -> int -> float -> float -> float array -> int -> int -> unit
(** [dlaset uplo m n alpha beta a ao lda] sets the strictly upper or strictly
    lower triangular or trapezoidal part of the m by n matrix [a], or the
    whole of it when [uplo] is [Full], to [alpha], and then sets the first
    [min m n] diagonal elements to [beta]. *)

val dlacpy : uplo -> int -> int -> float array -> int -> int -> float array -> int -> int -> unit
(** [dlacpy uplo m n a ao lda b bo ldb] copies the upper or lower triangular
    or trapezoidal part of the m by n matrix [a], or the whole of it when
    [uplo] is [Full], into [b]. *)

val dlange : norm -> int -> int -> float array -> int -> int -> float array -> float
(** [dlange norm m n a ao lda work] returns the selected norm of the m by n
    matrix [a], which is 0 when [m] or [n] is 0, and propagates a NaN
    element. [work] is a workspace of at least [m] elements and is read and
    written only when [norm] is [Inf_norm]. *)

val ieeeck : int -> float -> float -> int
(** [ieeeck ispec zero one] returns 1 when infinity arithmetic, and for
    [ispec] other than 0 also NaN arithmetic, behaves as IEEE requires
    without trapping, and 0 otherwise. [zero] and [one] are the values 0 and
    1, passed in so the checks cannot be folded away. *)

val iparmq : int -> string -> string -> int -> int -> int -> int -> int
(** [iparmq ispec name opts n ilo ihi lwork] returns a tuning parameter for
    the multishift QR algorithm, selected by [ispec] from 12 to 17, and -1
    for any other [ispec]. [name] is the calling routine, [ilo] and [ihi]
    bound the active block, and [n], [opts] and [lwork] are unused. *)

val ilaenv : int -> string -> string -> int -> int -> int -> int -> int
(** [ilaenv ispec name opts n1 n2 n3 n4] returns the blocking or tuning
    parameter selected by [ispec] from 1 to 17 for the routine [name], and -1
    for any other [ispec]. [name] is a Fortran style routine name such as
    "DGETRF", uppercased only when its first character is lower case, [opts]
    is the options string, and [n1] to [n4] are the problem dimensions. *)
