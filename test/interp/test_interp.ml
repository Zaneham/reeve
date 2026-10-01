(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Level 1, level 2 and level 3 checks for the polynomial interpolation
   routines, with the points, evaluation abscissas and tolerances taken from
   the slatec-modern Fortran suites, the level 3 numbers being the IBM
   System/360 values captured through Hercules. *)

open Reeve.Interp

let tol = 1.0e-10
let bad = ref 0

let fail what =
  Printf.printf "FAIL %s\n" what;
  bad := !bad + 1

let check what c = if not c then fail what
let near_t t what a b = if not (Float.abs (a -. b) < t) then fail what
let near what a b = near_t tol what a b

let same_t t what (a : float array) b =
  if Array.length a <> Array.length b then fail what
  else
    Array.iteri
      (fun i x -> near_t t (Printf.sprintf "%s[%d]" what i) x b.(i))
      a

let same what a b = same_t tol what a b
let ierr_is what (a : int) b = if a <> b then fail what

let table x y =
  let n = Array.length x in
  let c = Array.make n 0.0 in
  dplint n x y c;
  c

let value x c xx =
  let n = Array.length x in
  let yp = Array.make 1 0.0 in
  let work = Array.make (2 * n) 0.0 in
  let yfit, ierr = dpolvl 0 xx yp n x c work in
  ierr_is "dpolvl ierr" ierr 1;
  yfit

let derivs nder x c xx =
  let n = Array.length x in
  let yp = Array.make nder 0.0 in
  let work = Array.make (2 * n) 0.0 in
  let yfit, ierr = dpolvl nder xx yp n x c work in
  ierr_is "dpolvl ierr" ierr 1;
  (yfit, yp)

let coeffs xx x c =
  let n = Array.length x in
  let d = Array.make n 0.0 in
  let work = Array.make (2 * n) 0.0 in
  dpolcf xx n x c d work;
  d

let horner d xx z =
  let t = z -. xx in
  let s = ref 0.0 in
  for i = Array.length d - 1 downto 0 do
    s := (!s *. t) +. d.(i)
  done;
  !s

let raises what f =
  match f () with
  | () -> fail what
  | exception Invalid_argument _ -> ()

(* ---- Level 1, the regression suite ---- *)

let l1_dplint_quadratic () =
  let x = [| 0.0; 1.0; 2.0 |] and y = [| 0.0; 1.0; 4.0 |] in
  let c = table x y in
  near "l1 dplint quadratic at 1.5" (value x c 1.5) 2.25

let l1_dpolcf_taylor () =
  let x = [| -2.0; -1.0; 0.0; 1.0; 2.0 |] in
  let y = Array.map (fun v -> ((v *. v) -. 1.0) ** 2.0) x in
  let c = table x y in
  same_t 1.0e-8 "l1 dpolcf taylor" (coeffs 0.0 x c)
    [| 1.0; 0.0; -2.0; 0.0; 1.0 |]

let l1_dpolvl_derivatives () =
  let x = [| 0.0; 1.0; 2.0; 3.0 |] and y = [| 0.0; 1.0; 8.0; 27.0 |] in
  let c = table x y in
  let _, yp = derivs 3 x c 1.0 in
  same "l1 dpolvl derivatives of cube at one" yp [| 3.0; 6.0; 6.0 |]

(* ---- Level 2, the mathematical suite ---- *)

let l2_uniqueness () =
  let x1 = [| 0.0; 1.0; 2.0; 3.0; 4.0 |]
  and y1 = [| 0.0; 1.0; 8.0; 27.0; 64.0 |] in
  let x2 = [| 4.0; 2.0; 0.0; 3.0; 1.0 |]
  and y2 = [| 64.0; 8.0; 0.0; 27.0; 1.0 |] in
  let c1 = table x1 y1 and c2 = table x2 y2 in
  near "l2 uniqueness at 1.5" (value x1 c1 1.5) (value x2 c2 1.5);
  near "l2 uniqueness at 2.7" (value x1 c1 2.7) (value x2 c2 2.7);
  near "l2 uniqueness cube at 2.5" (value x1 c1 2.5) (2.5 ** 3.0);
  same_t 1.0e-8 "l2 uniqueness taylor agrees" (coeffs 0.0 x1 c1)
    (coeffs 0.0 x2 c2)

let l2_exactness f name =
  let x = Array.init 5 (fun i -> float_of_int i) in
  let y = Array.map f x in
  let c = table x y in
  for i = 1 to 20 do
    let xx = float_of_int (i - 1) /. 5.0 in
    near_t tol
      (Printf.sprintf "l2 exactness %s at %g" name xx)
      (value x c xx) (f xx)
  done

let l2_exactness_linear () = l2_exactness (fun v -> (2.0 *. v) +. 3.0) "linear"

let l2_exactness_quadratic () =
  l2_exactness (fun v -> (v -. 1.0) ** 2.0) "quadratic"

let l2_exactness_cubic () =
  l2_exactness
    (fun v -> (v ** 3.0) -. (3.0 *. (v ** 2.0)) +. (2.0 *. v))
    "cubic"

let l2_exactness_quartic () = l2_exactness (fun v -> v ** 4.0) "quartic"

(* ---- Level 3, the IBM System/360 values ---- *)

let l3_cube_ibm360 () =
  let t = 1.0e-12 in
  let x = Array.init 5 (fun i -> float_of_int i) in
  let y = Array.map (fun v -> v *. v *. v) x in
  let c = table x y in
  same_t t "l3 newton coefficients" c [| 0.0; 1.0; 3.0; 1.0; 0.0 |];
  near_t t "l3 p(0.5)" (value x c 0.5) 0.125;
  near_t t "l3 p(1.5)" (value x c 1.5) 3.375;
  near_t t "l3 p(2.5)" (value x c 2.5) 15.625;
  near_t t "l3 p(3.5)" (value x c 3.5) 42.875

(* ---- Reproduction of a polynomial of degree below the node count ---- *)

let p_cubic v = (2.0 *. v *. v *. v) -. (5.0 *. v *. v) +. v -. 7.0
let p_cubic_d1 v = (6.0 *. v *. v) -. (10.0 *. v) +. 1.0
let p_cubic_d2 v = (12.0 *. v) -. 10.0

let reproduce_cubic () =
  let x = [| -2.0; -1.0; 0.0; 1.0; 2.0; 3.0 |] in
  let y = Array.map p_cubic x in
  let c = table x y in
  for i = 0 to 40 do
    let xx = -2.0 +. (float_of_int i *. 0.125) in
    near (Printf.sprintf "reproduce cubic at %g" xx) (value x c xx) (p_cubic xx)
  done;
  let yfit, yp = derivs 5 x c 0.75 in
  near "reproduce cubic value at 0.75" yfit (p_cubic 0.75);
  same "reproduce cubic derivatives at 0.75" yp
    [| p_cubic_d1 0.75; p_cubic_d2 0.75; 12.0; 0.0; 0.0 |]

let fewer_derivatives_than_degree () =
  let x = [| -2.0; -1.0; 0.0; 1.0; 2.0; 3.0 |] in
  let y = Array.map p_cubic x in
  let c = table x y in
  let yfit, yp = derivs 2 x c (-0.5) in
  near "fewer derivatives value" yfit (p_cubic (-0.5));
  same "fewer derivatives" yp [| p_cubic_d1 (-0.5); p_cubic_d2 (-0.5) |];
  let yfit, yp = derivs 1 x c 2.25 in
  near "one derivative value" yfit (p_cubic 2.25);
  same "one derivative" yp [| p_cubic_d1 2.25 |]

let p_quartic v = (v ** 4.0) -. (3.0 *. (v ** 3.0)) +. (2.0 *. v) -. 1.0
let p_quartic_d1 v = (4.0 *. (v ** 3.0)) -. (9.0 *. (v ** 2.0)) +. 2.0
let p_quartic_d2 v = (12.0 *. (v ** 2.0)) -. (18.0 *. v)
let p_quartic_d3 v = (24.0 *. v) -. 18.0

let reproduce_quartic_uneven () =
  let x = [| 0.5; 1.25; 2.0; 3.5; 4.25 |] in
  let y = Array.map p_quartic x in
  let c = table x y in
  for i = 0 to 30 do
    let xx = 0.5 +. (float_of_int i *. 0.125) in
    near_t 1.0e-9
      (Printf.sprintf "reproduce quartic at %g" xx)
      (value x c xx) (p_quartic xx)
  done;
  let yfit, yp = derivs 4 x c 2.75 in
  near_t 1.0e-9 "reproduce quartic value at 2.75" yfit (p_quartic 2.75);
  same_t 1.0e-9 "reproduce quartic derivatives at 2.75" yp
    [| p_quartic_d1 2.75; p_quartic_d2 2.75; p_quartic_d3 2.75; 24.0 |]

let reproduce_nodes () =
  let x = [| -1.5; 0.25; 1.0; 2.75; 4.0; 5.5; 6.25 |] in
  let y = Array.map p_quartic x in
  let c = table x y in
  Array.iteri
    (fun i xi ->
      near_t 1.0e-9
        (Printf.sprintf "reproduce at node %d" i)
        (value x c xi) y.(i))
    x

(* ---- Taylor coefficients against the evaluation ---- *)

let taylor_round_trip () =
  let x = [| -1.0; 0.5; 1.5; 2.5; 4.0 |] in
  let y = Array.map (fun v -> Float.exp v) x in
  let c = table x y in
  List.iter
    (fun xx ->
      let d = coeffs xx x c in
      near_t 1.0e-9
        (Printf.sprintf "taylor value about %g" xx)
        d.(0) (value x c xx);
      List.iter
        (fun z ->
          near_t 1.0e-9
            (Printf.sprintf "taylor horner about %g at %g" xx z)
            (horner d xx z) (value x c z))
        [ -1.0; 0.0; 1.0; 2.0; 3.0; 4.0 ])
    [ 0.0; 1.0; 2.0 ]

let taylor_matches_derivatives () =
  let x = [| 0.0; 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let y = Array.map (fun v -> (v ** 5.0) -. (2.0 *. v) +. 3.0) x in
  let c = table x y in
  let xx = 1.5 in
  let d = coeffs xx x c in
  let _, yp = derivs 5 x c xx in
  let fac = ref 1.0 in
  for k = 1 to 5 do
    fac := !fac *. float_of_int k;
    near_t 1.0e-8
      (Printf.sprintf "taylor coefficient %d from derivative" k)
      (d.(k) *. !fac) yp.(k - 1)
  done

(* ---- Edge cases the Fortran singles out ---- *)

let one_point () =
  let x = [| 2.5 |] and y = [| 7.25 |] in
  let c = table x y in
  same "one point table" c [| 7.25 |];
  near "one point value" (value x c 100.0) 7.25;
  let yfit, yp = derivs 3 x c 100.0 in
  near "one point value with derivatives" yfit 7.25;
  same "one point derivatives" yp [| 0.0; 0.0; 0.0 |];
  same "one point taylor" (coeffs 3.0 x c) [| 7.25 |]

let two_points () =
  let x = [| 1.0; 3.0 |] and y = [| 5.0; 11.0 |] in
  let c = table x y in
  near "two points value" (value x c 2.0) 8.0;
  let yfit, yp = derivs 1 x c 2.0 in
  near "two points value with derivatives" yfit 8.0;
  same "two points slope" yp [| 3.0 |];
  let _, yp = derivs 4 x c 2.0 in
  same "two points excess derivatives" yp [| 3.0; 0.0; 0.0; 0.0 |];
  same "two points taylor about zero" (coeffs 0.0 x c) [| 2.0; 3.0 |]

let excess_derivatives () =
  let x = [| 0.0; 1.0; 2.0 |] and y = [| 1.0; 3.0; 9.0 |] in
  let c = table x y in
  let yfit, yp = derivs 5 x c 1.0 in
  near "excess derivatives value" yfit 3.0;
  same "excess derivatives" yp [| 4.0; 4.0; 0.0; 0.0; 0.0 |]

let work_untouched_at_nder_zero () =
  let x = [| 0.0; 1.0; 2.0 |] and y = [| 1.0; 3.0; 9.0 |] in
  let c = table x y in
  let work = Array.make 6 Float.nan in
  let yp = Array.make 1 Float.nan in
  let yfit, ierr = dpolvl 0 1.0 yp 3 x c work in
  ierr_is "work untouched ierr" ierr 1;
  near "work untouched value" yfit 3.0;
  check "work untouched" (Array.for_all Float.is_nan work);
  check "yp untouched" (Float.is_nan yp.(0))

let bad_arguments () =
  raises "dplint n zero" (fun () -> dplint 0 [| 1.0 |] [| 1.0 |] [| 0.0 |]);
  raises "dplint n negative" (fun () ->
      dplint (-1) [| 1.0 |] [| 1.0 |] [| 0.0 |]);
  raises "dplint repeated abscissa" (fun () ->
      dplint 3 [| 0.0; 1.0; 0.0 |] [| 1.0; 2.0; 3.0 |] (Array.make 3 0.0));
  raises "dplint repeated adjacent abscissa" (fun () ->
      dplint 2 [| 4.0; 4.0 |] [| 1.0; 2.0 |] (Array.make 2 0.0))

let () =
  l1_dplint_quadratic ();
  l1_dpolcf_taylor ();
  l1_dpolvl_derivatives ();
  l2_uniqueness ();
  l2_exactness_linear ();
  l2_exactness_quadratic ();
  l2_exactness_cubic ();
  l2_exactness_quartic ();
  l3_cube_ibm360 ();
  reproduce_cubic ();
  fewer_derivatives_than_degree ();
  reproduce_quartic_uneven ();
  reproduce_nodes ();
  taylor_round_trip ();
  taylor_matches_derivatives ();
  one_point ();
  two_points ();
  excess_derivatives ();
  work_untouched_at_nder_zero ();
  bad_arguments ();
  if !bad = 0 then print_string "interp: all checks passed\n" else exit 1
