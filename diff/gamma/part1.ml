(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Bessel functions, ported from the SLATEC double precision routines.
   Fixed order zero and one for I and K, and the Amos sequence routines
   for I, J and K of arbitrary non-negative order. *)

let tiny_dp = min_float
let huge_dp = max_float
let eps_2_dp = epsilon_float /. 2.0
let log10_radix_dp = log10 2.0
let digits_dp = 53
let min_exp_dp = -1021

let iter_limit = 1_000_000
