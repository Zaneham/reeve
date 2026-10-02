(* Reeve, Copyright 2026 Zane Hambly.
   The double precision data handling routines, ported from the SLATEC
   Fortran: the Singleton quicksort, the passive sort that returns a
   permutation, the routine that applies one, and the two way merge.
   Permutation vectors hold zero based indices, as pivot vectors do elsewhere
   in the library. Neither sort is stable: Singleton's quicksort moves equal
   keys past one another, so the order equal keys come out in is not the
   order they went in, and the permutations the passive sort returns are not
   the lexicographically least ones. The merge is stable, taking ties from the
   first run. *)

module Raw : sig
  (** The Fortran shapes, argument for argument. These are what the
      differential in [diff/] compares against the reference Fortran. Use the
      functions below unless you specifically want to pass [n] yourself. *)

  type kflag = Increasing | Increasing_carry | Decreasing | Decreasing_carry
  (** The Fortran's signed [KFLAG], with the four values it accepts as the four
      constructors. The sign chooses the order, so [Increasing] and
      [Increasing_carry] sort upwards and the other two downwards, and the
      magnitude chooses whether the second array moves: a [_carry] flag is the
      Fortran's 2 and -2, the plain one its 1 and -1. The Fortran's check that
      [KFLAG] is one of those four is gone, because nothing else can be
      written. *)

  val dsort : float array -> float array -> int -> kflag -> unit
  (** [dsort dx dy n kflag] sorts the first [n] elements of [dx] in place by
      Singleton's quicksort, making the same interchanges in [dy] when [kflag]
      is [Increasing_carry] or [Decreasing_carry]. For [Increasing] and
      [Decreasing], [dy] is never referenced and may be any array, including an
      empty one. Raises [Invalid_argument] when [n] is below one, where the
      Fortran reports through [XERMSG]. The explicit stack is 21 entries as it
      is in the Fortran, which is enough for [n] below 2^21 because the longer
      half of each split is the one stacked, and above that an index out of
      bounds escapes. *)

  val dpsort : float array -> int -> int array -> kflag -> unit
  (** [dpsort dx n iperm kflag] fills [iperm] with the permutation that sorts
      the first [n] elements of [dx], so that [dx.(iperm.(i))] is the [i]th
      value in sorted order, and leaves [dx] alone unless [kflag] is
      [Increasing_carry] or [Decreasing_carry], in which case [dx] is
      rearranged into that order as well. [iperm] holds zero based indices and
      is a valid permutation of [0] to [n - 1] on return either way. Raises
      [Invalid_argument] when [n] is below one. The Fortran's [IER] is gone:
      its two values both report a bad argument, and both now raise. The stack
      bound on [n] is the same as {!dsort}'s. *)

  val dpperm : float array -> int -> int array -> unit
  (** [dpperm dx n iperm] rearranges the first [n] elements of [dx] in place
      according to [iperm], so that element [i] afterwards is the element
      [iperm.(i)] from before. [iperm] holds zero based indices and is left
      unchanged. Raises [Invalid_argument] when [n] is below one, and when
      [iperm] is not a permutation of [0] to [n - 1], where the Fortran returns
      [IER] of 1 and 2. *)

  val d1merg : float array -> int -> int -> int -> int -> int -> unit
  (** [d1merg tcos i1 m1 i2 m2 i3] merges the ascending run of [m1] elements of
      [tcos] that starts at [i1] with the ascending run of [m2] elements that
      starts at [i2], writing the merged run of [m1 + m2] elements starting at
      [i3]. The three offsets are the Fortran's, which count elements skipped,
      and because the array is zero based here an offset is also the index of
      the element its run starts on. Ties go to the first run. The tail of
      whichever run is left over is copied forward element by element, as the
      Fortran's [DCOPY] does, so a destination overlapping a run above it will
      read what it has already written. [tcos] must be long enough for all
      three runs, which is the larger of [i1 + m1], [i2 + m2] and
      [i3 + m1 + m2]. *)
end

(* ---- The OCaml surface ---- *)

val sort : ?desc:bool -> ?carry:float array -> float array -> unit
(** [sort dx] sorts the whole of [dx] upwards in place, or downwards with
    [~desc:true]. [~carry:dy] makes the same interchanges in [dy], which has to
    be at least as long as [dx] or you get [Invalid_argument]. The sort isn't
    stable, so equal keys come out in whatever order the partitioning left
    them. *)

val sort_index : ?desc:bool -> float array -> int array
(** [sort_index dx] returns the permutation that sorts [dx] without touching
    [dx] itself, so [dx.(i)] for [i] taken from the result comes out in order.
    [~desc:true] sorts downwards. It isn't the lexicographically least
    permutation when there are ties. *)

val permute : float array -> int array -> unit
(** [permute dx iperm] rearranges [dx] in place so element [i] afterwards is
    the element [iperm.(i)] from before, which is what you want to replay a
    {!sort_index} onto another array. [iperm] has to be a permutation of [0] to
    [Array.length dx - 1] and is left alone. *)

val merge : float array -> float array -> float array
(** [merge a b] merges two arrays that are already ascending into a fresh one
    of their combined length. It's stable, so a tie takes the element from [a]
    first. Nothing checks that the inputs are sorted. *)
