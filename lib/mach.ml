(* Reeve, Copyright 2026 Zane Hambly.
   The machine constants the ported routines ask D1MACH and DLAMCH for, the
   float32 rounding the Fortran gets from a single precision expression, and
   the Blue scaling thresholds that follow from the exponent range. This file
   depends on nothing. *)

(* ---- D1MACH ---- *)

let eps_dp = epsilon_float
let eps_2_dp = epsilon_float /. 2.0
let tiny_dp = min_float
let huge_dp = max_float
let log10_radix_dp = log10 2.0
let digits_dp = 53
let min_exp_dp = -1021
let max_exp_dp = 1024

(* ---- Single precision rounding ---- *)

let single x = Int32.float_of_bits (Int32.bits_of_float x)

(* ---- DLAMCH ---- *)

let fdigits = 2 - snd (Float.frexp epsilon_float)
let fminexp = snd (Float.frexp min_float)
let fmaxexp = snd (Float.frexp max_float)

(* ---- Blue's scaling thresholds ---- *)

let tsml =
  Float.ldexp 1.0 (int_of_float (Float.ceil (float_of_int (fminexp - 1) *. 0.5)))

let tbig =
  Float.ldexp 1.0
    (int_of_float (Float.floor (float_of_int (fmaxexp - fdigits + 1) *. 0.5)))

let ssml =
  Float.ldexp 1.0
    (-int_of_float (Float.floor (float_of_int (fminexp - fdigits) *. 0.5)))

let sbig =
  Float.ldexp 1.0
    (-int_of_float (Float.ceil (float_of_int (fmaxexp + fdigits - 1) *. 0.5)))
