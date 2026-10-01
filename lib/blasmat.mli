(* Reeve, Copyright 2026 Zane Hambly.
   BLAS level 2 and level 3 for dense matrices, ported from the reference
   Fortran. A matrix is a [float array] plus a base offset [ao] and a leading
   dimension [lda], column major, so the element (i, j) zero based lives at
   [ao + i + j * lda]. The Fortran character options are variants here. *)

type trans = No_trans | Trans | Conj_trans
(** The [TRANS] option of the Fortran. [Conj_trans] is [Trans] for real data,
    as it is in the Fortran. *)

type uplo = Upper | Lower
(** The [UPLO] option, saying which triangle of a matrix is referenced. *)

type side = Left | Right
(** The [SIDE] option, saying whether the triangular matrix multiplies from
    the left or the right. *)

type diag = Unit | Non_unit
(** The [DIAG] option. [Unit] means the diagonal of the triangular matrix is
    taken as one and not referenced. *)

val dgemv : trans -> int -> int -> float -> float array -> int -> int
  -> float array -> int -> int -> float -> float array -> int -> int -> unit
(** [dgemv trans m n alpha a ao lda x xo incx beta y yo incy] forms
    [y <- alpha * op a * x + beta * y], where [op a] is the [m] by [n] matrix
    [a] when [trans] is [No_trans] and its transpose otherwise. [x] has [n]
    elements and [y] has [m] when [trans] is [No_trans], the other way round
    otherwise. Raises [Invalid_argument] naming the routine and the Fortran
    argument position when [m] or [n] is negative, [lda] is below
    [max 1 m], or either increment is zero. *)

val dger : int -> int -> float -> float array -> int -> int -> float array
  -> int -> int -> float array -> int -> int -> unit
(** [dger m n alpha x xo incx y yo incy a ao lda] adds the rank one update
    [alpha * x * transpose y] to the [m] by [n] matrix [a]. Raises
    [Invalid_argument] naming the routine and the Fortran argument position
    when [m] or [n] is negative, either increment is zero, or [lda] is below
    [max 1 m]. *)

val dgemm : trans -> trans -> int -> int -> int -> float -> float array -> int
  -> int -> float array -> int -> int -> float -> float array -> int -> int
  -> unit
(** [dgemm transa transb m n k alpha a ao lda b bo ldb beta c co ldc] forms
    [c <- alpha * op a * op b + beta * c], where [c] is [m] by [n], [op a] is
    [m] by [k] and [op b] is [k] by [n], each [op] being the transpose when
    its [trans] is not [No_trans]. Raises [Invalid_argument] naming the
    routine and the Fortran argument position when [m], [n] or [k] is
    negative or a leading dimension is too small for the shape it carries. *)

val dsyrk : uplo -> trans -> int -> int -> float -> float array -> int -> int
  -> float -> float array -> int -> int -> unit
(** [dsyrk uplo trans n k alpha a ao lda beta c co ldc] forms the symmetric
    rank [k] update [c <- alpha * a * transpose a + beta * c] when [trans] is
    [No_trans], with [a] of shape [n] by [k], and
    [c <- alpha * transpose a * a + beta * c] otherwise, with [a] of shape [k]
    by [n]. Only the [uplo] triangle of the order [n] matrix [c] is read and
    written. Raises [Invalid_argument] naming the routine and the Fortran
    argument position when [n] or [k] is negative or a leading dimension is
    too small. *)

val dtrsm : side -> uplo -> trans -> diag -> int -> int -> float
  -> float array -> int -> int -> float array -> int -> int -> unit
(** [dtrsm side uplo transa diag m n alpha a ao lda b bo ldb] solves in place
    in the [m] by [n] matrix [b] the triangular system
    [op a * x = alpha * b] when [side] is [Left] and
    [x * op a = alpha * b] when it is [Right], where [a] is triangular of
    order [m] or [n] respectively, [uplo] says which triangle holds it,
    [diag] says whether its diagonal is taken as one, and [op a] is its
    transpose unless [transa] is [No_trans]. Raises [Invalid_argument] naming
    the routine and the Fortran argument position when [m] or [n] is negative
    or a leading dimension is too small. *)

val dtrmm : side -> uplo -> trans -> diag -> int -> int -> float
  -> float array -> int -> int -> float array -> int -> int -> unit
(** [dtrmm side uplo transa diag m n alpha a ao lda b bo ldb] overwrites the
    [m] by [n] matrix [b] with [alpha * op a * b] when [side] is [Left] and
    [alpha * b * op a] when it is [Right], the triangular matrix [a] and the
    options reading as in {!dtrsm}. Raises [Invalid_argument] naming the
    routine and the Fortran argument position when [m] or [n] is negative or
    a leading dimension is too small. *)
