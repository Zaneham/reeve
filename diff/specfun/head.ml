(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Gamma, psi, the Airy functions and the Carlson elliptic integrals, ported
   from the SLATEC double precision routines. *)

let eps_2_dp = epsilon_float /. 2.0
let eps_dp = epsilon_float
let tiny_dp = min_float
let huge_dp = max_float
let log10_radix_dp = log10 2.0

let single x = Int32.float_of_bits (Int32.bits_of_float x)
let cube x = x *. x *. x

let csevl x cs n =
  if Float.abs x > 1.0 then
    invalid_arg "csevl: x outside the interval (-1,+1)";
  let twox = 2.0 *. x in
  let b0 = ref 0.0 and b1 = ref 0.0 and b2 = ref 0.0 in
  for i = 1 to n do
    b2 := !b1;
    b1 := !b0;
    b0 := (twox *. !b1) -. !b2 +. cs.(n - i)
  done;
  0.5 *. (!b0 -. !b2)

