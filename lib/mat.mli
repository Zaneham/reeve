(* Reeve, Copyright 2026 Zane Hambly.
   The dense matrix the library hands around. *)

type t = { data : float array; rows : int; cols : int }
(** Column major, so element [(i, j)] zero based is [data.(i + j * rows)].
    The record is open so you can hand [data] to something else, but prefer
    {!get} and {!set}. *)

val create : int -> int -> t
(** [create rows cols] is zeroed. *)

val of_array : int -> int -> float array -> t
(** [of_array rows cols data] takes [data] as the storage, without copying. *)

val of_vec : float array -> t
(** [of_vec v] is [v] as a single column, sharing the storage. *)

val init : int -> int -> (int -> int -> float) -> t
(** [init rows cols f] fills from [f i j]. *)

val copy : t -> t

val rows : t -> int
val cols : t -> int
val data : t -> float array

val get : t -> int -> int -> float
val set : t -> int -> int -> float -> unit

val square : string -> t -> unit
(** [square name m] raises [Invalid_argument] naming [name] unless [m] is
    square. For the routines that need it. *)

val same_rows : string -> t -> t -> unit
(** [same_rows name a b] raises [Invalid_argument] naming [name] unless [a]
    and [b] have the same number of rows. *)
