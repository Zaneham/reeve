(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Level 1, level 2 and level 3 checks for gamma, psi, the Airy functions and
   the Carlson elliptic integrals. *)

open Reeve.Specfun

let pi = 3.14159265358979323846
let euler_gamma = 0.5772156649015328606
let tol = 1.0e-10
let loose_tol = 1.0e-6
let hist_tol = 1.0e-8
let bad = ref 0

let fail what =
  Printf.printf "FAIL %s\n" what;
  bad := !bad + 1

let near what a b = if not (Float.abs (a -. b) < tol) then fail what
let loose what a b = if not (Float.abs (a -. b) < loose_tol) then fail what
let near_t what t a b = if not (Float.abs (a -. b) < t) then fail what

let rel_t what t a b =
  let e =
    if Float.abs b > 1.0e-14 then Float.abs (a -. b) /. Float.abs b
    else Float.abs (a -. b)
  in
  if not (e < t) then fail what

let holds what c = if not c then fail what

(* ---- Level 1, regression ---- *)

let l1_gamma () =
  near "l1 dgamln 1" (dgamln 1.0) 0.0;
  near "l1 dgamln 2" (dgamln 2.0) 0.0;
  near "l1 dgamln 0.5" (dgamln 0.5) (0.5 *. log pi);
  loose "l1 dpsi 1" (dpsi 1.0) (-.euler_gamma)

let l1_airy () =
  loose "l1 dai 0" (dai 0.0) 0.35502805388781724;
  loose "l1 dbi 0" (dbi 0.0) 0.61492662744600073;
  loose "l1 dai -1" (dai (-1.0)) 0.5355608832923521;
  holds "l1 dai decays for x > 0"
    (dai 1.0 > 0.0 && dai 5.0 > 0.0 && dai 5.0 < dai 1.0);
  holds "l1 dbi grows for x > 0" (dbi 1.0 > 0.0 && dbi 5.0 > dbi 1.0)

let l1_elliptic () =
  near "l1 drf 1 1 1" (drf 1.0 1.0 1.0) 1.0;
  holds "l1 drf 0 1 2 in 1.3 to 1.4"
    (drf 0.0 1.0 2.0 > 1.3 && drf 0.0 1.0 2.0 < 1.4);
  near "l1 drc 1 1" (drc 1.0 1.0) 1.0;
  loose "l1 drc 0 1" (drc 0.0 1.0) (pi /. 2.0);
  near "l1 drd 1 1 1" (drd 1.0 1.0 1.0) 1.0

(* ---- Level 2, mathematical identities ---- *)

let l2_gamma () =
  let x = 3.5 in
  near "l2 dgamln recurrence" (dgamln (x +. 1.0)) (log x +. dgamln x);
  let x = 0.3 in
  loose "l2 dgamln reflection"
    (dgamln x +. dgamln (1.0 -. x))
    (log pi -. log (Float.abs (sin (pi *. x))));
  let x = 1.5 in
  loose "l2 dgamln duplication"
    (dgamln (2.0 *. x))
    (((2.0 *. x) -. 1.0) *. log 2.0 +. dgamln x +. dgamln (x +. 0.5)
     -. (0.5 *. log pi));
  let x = 2.5 in
  loose "l2 dpsi recurrence" (dpsi (x +. 1.0)) (dpsi x +. (1.0 /. x));
  let x = 0.25 in
  near_t "l2 dpsi reflection" 1.0e-3
    (dpsi (1.0 -. x) -. dpsi x)
    (pi *. cos (pi *. x) /. sin (pi *. x))

let l2_airy () =
  let h = 1.0e-6 in
  let aip = (dai (0.0 +. h) -. dai (0.0 -. h)) /. (2.0 *. h) in
  let bip = (dbi (0.0 +. h) -. dbi (0.0 -. h)) /. (2.0 *. h) in
  near_t "l2 wronskian ai bi" 0.01
    ((dai 0.0 *. bip) -. (aip *. dbi 0.0))
    (1.0 /. pi);
  near_t "l2 dbi 0 over dai 0" 0.01 (dbi 0.0 /. dai 0.0) (sqrt 3.0);
  let x = 1.0 and h = 1.0e-4 in
  let second =
    (dai (x +. h) -. (2.0 *. dai x) +. dai (x -. h)) /. (h *. h)
  in
  near_t "l2 dai ode" 0.01 second (x *. dai x);
  let x = 5.0 in
  let asymp =
    exp (-2.0 /. 3.0 *. (x ** 1.5)) /. (2.0 *. sqrt pi *. (x ** 0.25))
  in
  near_t "l2 dai asymptotic" 0.1 (dai x /. asymp) 1.0

let l2_elliptic () =
  let rf1 = drf 1.0 2.0 3.0 in
  let rf2 = drf 2.0 1.0 3.0 in
  let rf3 = drf 3.0 2.0 1.0 in
  near "l2 drf symmetry xy" rf1 rf2;
  near "l2 drf symmetry xz" rf2 rf3;
  let k = 4.0 in
  loose "l2 drf homogeneity" (drf (k *. 1.0) (k *. 2.0) (k *. 3.0) *. sqrt k) rf1;
  let x = 4.0 in
  near "l2 drc x x" (drc x x) (1.0 /. sqrt x);
  let rd1 = drd 1.0 2.0 3.0 in
  loose "l2 drd homogeneity"
    (drd (k *. 1.0) (k *. 2.0) (k *. 3.0) *. (k ** 1.5))
    rd1;
  loose "l2 drf 0 1 1" (drf 0.0 1.0 1.0) (pi /. 2.0)

(* ---- Level 3, classical tabulated values ---- *)

let l3_gamma () =
  rel_t "l3 dgamln 1.5" hist_tol (dgamln 1.5) (-0.12078223763524522);
  rel_t "l3 dgamln 2.5" hist_tol (dgamln 2.5) 0.28468287047291918;
  rel_t "l3 dgamln 5" tol (dgamln 5.0) 3.1780538303479458;
  rel_t "l3 dpsi 1" hist_tol (dpsi 1.0) (-0.5772156649015329);
  rel_t "l3 dpsi 2" hist_tol (dpsi 2.0) 0.4227843350984671

let l3_airy () =
  rel_t "l3 dai 0" hist_tol (dai 0.0) 0.35502805388781724;
  rel_t "l3 dbi 0" hist_tol (dbi 0.0) 0.61492662744600073;
  rel_t "l3 dai -1" hist_tol (dai (-1.0)) 0.5355608832923521;
  rel_t "l3 dai 1" 0.01 (dai 1.0) 0.1352924163128814;
  rel_t "l3 dbi 1" 0.01 (dbi 1.0) 1.2074235949528713

let l3_elliptic () =
  rel_t "l3 drf 0 1 1" tol (drf 0.0 1.0 1.0) 1.5707963267948966;
  rel_t "l3 drf 1 1 1" tol (drf 1.0 1.0 1.0) 1.0;
  rel_t "l3 drc 1 2" hist_tol (drc 1.0 2.0) 0.7853981633974483;
  rel_t "l3 drf 0 0.5 1" hist_tol (drf 0.0 0.5 1.0) 1.8540746773013719;
  rel_t "l3 drd 0 2 1" 0.001 (drd 0.0 2.0 1.0) 1.7972103521033882

(* ---- Domains the Fortran refuses ---- *)

let raises what f =
  match f () with
  | _ -> fail what
  | exception Invalid_argument _ -> ()

let l1_domains () =
  raises "dgamln 0 raises" (fun () -> dgamln 0.0);
  raises "dgamln -1 raises" (fun () -> dgamln (-1.0));
  raises "dpsi 0 raises" (fun () -> dpsi 0.0);
  raises "dpsi -1 raises" (fun () -> dpsi (-1.0));
  raises "dpsi -2 raises" (fun () -> dpsi (-2.0));
  raises "dbi overflow raises" (fun () -> dbi 1.0e6);
  raises "drf negative raises" (fun () -> drf (-1.0) 1.0 1.0);
  raises "drd negative raises" (fun () -> drd (-1.0) 1.0 1.0);
  raises "drc negative raises" (fun () -> drc (-1.0) 1.0);
  raises "drc y zero raises" (fun () -> drc 1.0 0.0);
  raises "drf uplim raises" (fun () -> drf 0.0 1.0 Float.max_float);
  raises "drd uplim raises" (fun () -> drd 0.0 1.0 Float.max_float)

(* ---- The Fortran shapes ---- *)

let raw_dcsevl () =
  let cs = [| 0.5; 0.25; -0.125; 0.0625 |] in
  let t x k = Raw.dcsevl x cs k in
  List.iter
    (fun x ->
      let t0 = 1.0 and t1 = x in
      let t2 = (2.0 *. x *. x) -. 1.0 in
      let t3 = (4.0 *. x *. x *. x) -. (3.0 *. x) in
      near_t "raw dcsevl one term" tol (t x 1) (0.5 *. cs.(0) *. t0);
      near_t "raw dcsevl two terms" tol (t x 2)
        ((0.5 *. cs.(0) *. t0) +. (cs.(1) *. t1));
      near_t "raw dcsevl three terms" tol (t x 3)
        ((0.5 *. cs.(0) *. t0) +. (cs.(1) *. t1) +. (cs.(2) *. t2));
      near_t "raw dcsevl four terms" tol (t x 4)
        ((0.5 *. cs.(0) *. t0) +. (cs.(1) *. t1) +. (cs.(2) *. t2)
        +. (cs.(3) *. t3)))
    [ -1.0; -0.75; -0.25; 0.0; 0.3; 0.8; 1.0 ];
  near "raw dcsevl zero terms" (t 0.5 0) 0.0;
  raises "raw dcsevl above one raises" (fun () -> ignore (t 1.5 2));
  raises "raw dcsevl below minus one raises" (fun () -> ignore (t (-1.5) 2))

let raw_dcot () =
  List.iter
    (fun x ->
      rel_t "raw dcot against the libm" 1.0e-13 (Raw.dcot x) (1.0 /. tan x);
      holds "raw dcot is odd" (Raw.dcot (-.x) = -.Raw.dcot x))
    [ 0.125; 0.5; 1.0; 1.5; 2.0; 3.0; 4.0; 7.0 ];
  rel_t "raw dcot at a large argument" 1.0e-9 (Raw.dcot 1000.0)
    (1.0 /. tan 1000.0);
  rel_t "raw dcot has period pi" 1.0e-12 (Raw.dcot (0.7 +. pi)) (Raw.dcot 0.7);
  rel_t "raw dcot double angle" 1.0e-12 (Raw.dcot 1.4)
    ((Raw.dcot 0.7 -. (1.0 /. Raw.dcot 0.7)) /. 2.0);
  raises "raw dcot zero raises" (fun () -> ignore (Raw.dcot 0.0));
  raises "raw dcot huge argument raises" (fun () ->
      ignore (Raw.dcot 1.0e300))

let () =
  l1_gamma ();
  l1_airy ();
  l1_elliptic ();
  l1_domains ();
  l2_gamma ();
  l2_airy ();
  l2_elliptic ();
  l3_gamma ();
  l3_airy ();
  l3_elliptic ();
  raw_dcsevl ();
  raw_dcot ();
  if !bad = 0 then print_string "specfun: all checks passed\n" else exit 1
