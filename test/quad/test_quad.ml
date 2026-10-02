(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Checks for the quadrature routines. The DGAUS8 integrands, tolerances and
   expected values are the ones in the slatec-modern level 1 and level 2
   diff_integ suites; the rest are closed forms, or polynomials inside the
   degree each rule integrates exactly, so no reference number is needed. *)

open Reeve.Quad.Raw
module Quad = Reeve.Quad

let bad = ref 0

let fail what =
  Printf.printf "FAIL %s\n" what;
  bad := !bad + 1

let check what c = if not c then fail what
let near_t t what a b = if not (Float.abs (a -. b) < t) then fail what
let ierr_is what (a : int) b = if a <> b then fail what
let pi = 3.14159265358979323846

(* ---- DGAUS8, the slatec-modern level 1 and level 2 cases ---- *)

let l1_gaus8_x_squared () =
  let ans, ierr, _ = dgaus8 (fun x -> x *. x) 0.0 1.0 1.0e-12 in
  near_t 1.0e-10 "l1 dgaus8 x^2 over 0 1" ans (1.0 /. 3.0);
  ierr_is "l1 dgaus8 x^2 ierr" ierr 1

let l1_gaus8_sin () =
  let ans, ierr, _ = dgaus8 sin 0.0 pi 1.0e-12 in
  near_t 1.0e-10 "l1 dgaus8 sin over 0 pi" ans 2.0;
  ierr_is "l1 dgaus8 sin ierr" ierr 1

let l1_gaus8_exp_neg () =
  let ans, ierr, _ = dgaus8 (fun x -> exp (-.x)) 0.0 1.0 1.0e-12 in
  near_t 1.0e-10 "l1 dgaus8 exp(-x) over 0 1" ans (1.0 -. exp (-1.0));
  ierr_is "l1 dgaus8 exp(-x) ierr" ierr 1

let l2_gaus8_lorentzian () =
  let ans, ierr, _ =
    dgaus8 (fun x -> 1.0 /. (1.0 +. (x *. x))) 0.0 1.0 1.0e-12
  in
  near_t 1.0e-12 "l2 dgaus8 1/(1+x^2) over 0 1" ans (pi /. 4.0);
  ierr_is "l2 dgaus8 1/(1+x^2) ierr" ierr 1

(* ---- DGAUS8 against the rule and the status codes ---- *)

let gaus8_degree_fifteen () =
  let ans, ierr, _ = dgaus8 (fun x -> x ** 14.0) (-1.0) 1.0 1.0e-12 in
  near_t 1.0e-15 "dgaus8 x^14 over -1 1" ans (2.0 /. 15.0);
  ierr_is "dgaus8 x^14 ierr" ierr 1;
  let ans, ierr, _ = dgaus8 (fun x -> x ** 15.0) (-1.0) 1.0 1.0e-12 in
  near_t 1.0e-15 "dgaus8 x^15 over -1 1" ans 0.0;
  ierr_is "dgaus8 x^15 ierr" ierr 1

let gaus8_reversed () =
  let ans, ierr, _ = dgaus8 sin pi 0.0 1.0e-12 in
  near_t 1.0e-10 "dgaus8 sin over pi 0" ans (-2.0);
  ierr_is "dgaus8 reversed ierr" ierr 1

let gaus8_equal_limits () =
  let ans, ierr, e = dgaus8 sin 1.0 1.0 1.0e-12 in
  check "dgaus8 a = b ans is zero" (ans = 0.0);
  ierr_is "dgaus8 a = b ierr" ierr 1;
  check "dgaus8 a = b leaves a non-negative err alone" (e = 1.0e-12)

let gaus8_too_nearly_equal () =
  let ans, ierr, e = dgaus8 sin (1.0 -. 1.0e-15) 1.0 (-1.0) in
  check "dgaus8 nearly equal ans is zero" (ans = 0.0);
  ierr_is "dgaus8 nearly equal ierr" ierr (-1);
  check "dgaus8 nearly equal err estimate is zero" (e = 0.0)

let gaus8_error_estimate () =
  let cases =
    [ ("exp over 0 2", exp, 0.0, 2.0, exp 2.0 -. 1.0);
      ("sin over 0 pi", sin, 0.0, pi, 2.0);
      ("1/x over 1 10", (fun x -> 1.0 /. x), 1.0, 10.0, log 10.0);
      ("sqrt over 0 1", sqrt, 0.0, 1.0, 2.0 /. 3.0) ]
  in
  List.iter
    (fun (what, f, a, b, exact) ->
      let ans, _, e = dgaus8 f a b (-1.0e-10) in
      check
        (Printf.sprintf "dgaus8 error estimate bounds %s" what)
        (Float.abs (ans -. exact) <= Float.abs e +. 1.0e-13))
    cases

let gaus8_tighter_is_better () =
  let exact = log 10.0 in
  let coarse, _, _ = dgaus8 (fun x -> 1.0 /. x) 1.0 10.0 1.0e-4 in
  let fine, _, _ = dgaus8 (fun x -> 1.0 /. x) 1.0 10.0 1.0e-14 in
  check "dgaus8 a tighter tolerance is no worse"
    (Float.abs (fine -. exact) <= Float.abs (coarse -. exact) +. 1.0e-15);
  near_t 1.0e-13 "dgaus8 1/x over 1 10" fine exact

let inv_sqrt x = if x > 0.0 then 1.0 /. sqrt x else 0.0

let gaus8_insufficient_accuracy () =
  let ans, ierr, _ = dgaus8 inv_sqrt 0.0 1.0 1.0e-14 in
  ierr_is "dgaus8 reports it missed the tolerance" ierr 2;
  near_t 1.0e-4 "dgaus8 1/sqrt(x) over 0 1 is still close" ans 2.0

(* ---- DQNC79 ---- *)

let qnc79_degree_seven () =
  let ans, ierr, _ = dqnc79 (fun x -> x *. x *. x) 0.0 1.0 1.0e-12 in
  near_t 1.0e-15 "dqnc79 x^3 over 0 1" ans 0.25;
  ierr_is "dqnc79 x^3 ierr" ierr 1;
  let ans, ierr, _ = dqnc79 (fun x -> x ** 7.0) 0.0 1.0 1.0e-12 in
  near_t 1.0e-15 "dqnc79 x^7 over 0 1" ans 0.125;
  ierr_is "dqnc79 x^7 ierr" ierr 1

let qnc79_closed_forms () =
  let ans, ierr, _ = dqnc79 sin 0.0 pi 1.0e-12 in
  near_t 1.0e-11 "dqnc79 sin over 0 pi" ans 2.0;
  ierr_is "dqnc79 sin ierr" ierr 1;
  let ans, ierr, _ = dqnc79 (fun x -> exp (-.x)) 0.0 1.0 1.0e-12 in
  near_t 1.0e-12 "dqnc79 exp(-x) over 0 1" ans (1.0 -. exp (-1.0));
  ierr_is "dqnc79 exp(-x) ierr" ierr 1;
  let ans, ierr, _ =
    dqnc79 (fun x -> 1.0 /. (1.0 +. (x *. x))) 0.0 1.0 1.0e-12
  in
  near_t 1.0e-12 "dqnc79 1/(1+x^2) over 0 1" ans (pi /. 4.0);
  ierr_is "dqnc79 1/(1+x^2) ierr" ierr 1

let qnc79_counts_evaluations () =
  let seen = ref 0 in
  let f x =
    incr seen;
    exp (-.(x *. x))
  in
  let _, ierr, k = dqnc79 f 0.0 2.0 1.0e-12 in
  ierr_is "dqnc79 count ierr" ierr 1;
  check "dqnc79 reports at least thirteen evaluations" (k >= 13);
  check "dqnc79 count matches the calls made" (k = !seen)

let qnc79_reversed () =
  let ans, ierr, _ = dqnc79 sin pi 0.0 1.0e-12 in
  near_t 1.0e-11 "dqnc79 sin over pi 0" ans (-2.0);
  ierr_is "dqnc79 reversed ierr" ierr 1

let qnc79_equal_limits () =
  let ans, ierr, k = dqnc79 sin 1.0 1.0 1.0e-12 in
  check "dqnc79 a = b ans is zero" (ans = 0.0);
  ierr_is "dqnc79 a = b ierr" ierr (-1);
  check "dqnc79 a = b made no evaluations" (k = 0)

let qnc79_too_nearly_equal () =
  let ans, ierr, _ = dqnc79 sin (1.0 -. 1.0e-15) 1.0 1.0e-12 in
  check "dqnc79 nearly equal ans is zero" (ans = 0.0);
  ierr_is "dqnc79 nearly equal ierr" ierr (-1)

let qnc79_additive () =
  let f x = exp (-.x) *. cos x in
  let whole, _, _ = dqnc79 f 0.0 3.0 1.0e-13 in
  let left, _, _ = dqnc79 f 0.0 1.25 1.0e-13 in
  let right, _, _ = dqnc79 f 1.25 3.0 1.0e-13 in
  near_t 1.0e-12 "dqnc79 splits additively" whole (left +. right)

let qnc79_insufficient_accuracy () =
  let ans, ierr, k = dqnc79 inv_sqrt 0.0 1.0 1.0e-14 in
  ierr_is "dqnc79 reports it missed the tolerance" ierr 2;
  near_t 1.0e-5 "dqnc79 1/sqrt(x) over 0 1 is still close" ans 2.0;
  check "dqnc79 called the integrand many times on a hard problem" (k > 1000)

(* ---- DAVINT ---- *)

let quadratic x = 1.0 +. (2.0 *. x) +. (3.0 *. x *. x)
let quadratic_integral x = x +. (x *. x) +. (x *. x *. x)

let davint_parabola_exact () =
  let x = Array.init 5 (fun i -> float_of_int i *. 0.25) in
  let y = Array.map quadratic x in
  let ans, ierr = davint x y 5 0.0 1.0 in
  ierr_is "davint parabola ierr" ierr 1;
  near_t 1.0e-14 "davint 1+2x+3x^2 over 0 1" ans
    (quadratic_integral 1.0 -. quadratic_integral 0.0)

let davint_interior_limits () =
  let x = Array.init 9 (fun i -> float_of_int i *. 0.125) in
  let y = Array.map quadratic x in
  let ans, ierr = davint x y 9 0.1 0.9 in
  ierr_is "davint interior limits ierr" ierr 1;
  near_t 1.0e-13 "davint 1+2x+3x^2 over 0.1 0.9" ans
    (quadratic_integral 0.9 -. quadratic_integral 0.1)

let davint_uneven_abscissas () =
  let x = [| 0.0; 0.1; 0.35; 0.4; 0.72; 0.9; 1.0 |] in
  let y = Array.map quadratic x in
  let ans, ierr = davint x y 7 0.0 1.0 in
  ierr_is "davint uneven ierr" ierr 1;
  near_t 1.0e-13 "davint uneven abscissas" ans
    (quadratic_integral 1.0 -. quadratic_integral 0.0)

let davint_trapezoid () =
  let x = [| 0.0; 1.0 |] and y = [| 3.0; 5.0 |] in
  let ans, ierr = davint x y 2 0.0 1.0 in
  ierr_is "davint trapezoid ierr" ierr 1;
  near_t 1.0e-14 "davint trapezoid over the table" ans 4.0;
  let ans, ierr = davint x y 2 0.5 2.0 in
  ierr_is "davint trapezoid extrapolated ierr" ierr 1;
  near_t 1.0e-14 "davint trapezoid past the table" ans 8.25

let davint_equal_limits () =
  let x = Array.init 5 (fun i -> float_of_int i *. 0.25) in
  let y = Array.map quadratic x in
  let ans, ierr = davint x y 5 0.5 0.5 in
  check "davint xlo = xup ans is zero" (ans = 0.0);
  ierr_is "davint xlo = xup ierr" ierr 1

let davint_bad_bounds () =
  let x = Array.init 5 (fun i -> float_of_int i *. 0.25) in
  let y = Array.map quadratic x in
  let ans, ierr = davint x y 5 0.75 0.25 in
  check "davint xup below xlo ans is zero" (ans = 0.0);
  ierr_is "davint xup below xlo ierr" ierr 2

let davint_too_few_between () =
  let x = [| 0.0; 1.0; 2.0; 3.0 |] in
  let y = Array.map quadratic x in
  let ans, ierr = davint x y 4 2.5 3.0 in
  check "davint few points high ans is zero" (ans = 0.0);
  ierr_is "davint few points high ierr" ierr 3;
  let ans, ierr = davint x y 4 0.0 0.5 in
  check "davint few points low ans is zero" (ans = 0.0);
  ierr_is "davint few points low ierr" ierr 3

let davint_not_increasing () =
  let x = [| 0.0; 1.0; 1.0; 3.0 |] in
  let y = Array.map quadratic x in
  let ans, ierr = davint x y 4 0.0 3.0 in
  check "davint repeated abscissa ans is zero" (ans = 0.0);
  ierr_is "davint repeated abscissa ierr" ierr 4;
  let x = [| 0.0; 2.0; 1.0; 3.0 |] in
  let y = Array.map quadratic x in
  let ans, ierr = davint x y 4 0.0 3.0 in
  check "davint decreasing abscissa ans is zero" (ans = 0.0);
  ierr_is "davint decreasing abscissa ierr" ierr 4

let davint_too_few_values () =
  let x = [| 0.0 |] and y = [| 1.0 |] in
  let ans, ierr = davint x y 1 0.0 1.0 in
  check "davint n = 1 ans is zero" (ans = 0.0);
  ierr_is "davint n = 1 ierr" ierr 5

let davint_agrees_with_dgaus8 () =
  let f x = exp (-.x) *. cos (3.0 *. x) in
  let n = 401 in
  let x = Array.init n (fun i -> float_of_int i *. 2.0 /. float_of_int (n - 1)) in
  let y = Array.map f x in
  let ans, ierr = davint x y n 0.0 2.0 in
  ierr_is "davint against dgaus8 ierr" ierr 1;
  let want, _, _ = dgaus8 f 0.0 2.0 1.0e-13 in
  near_t 1.0e-8 "davint against dgaus8" ans want

(* ---- DPPGQ8 ---- *)

let pp_one_piece = [| 1.0; 2.0; 6.0 |]
let pp_breaks = [| 0.0; 1.0 |]

let ppgq8_spline_value () =
  let ans, ierr, _ =
    dppgq8 (fun _ -> 1.0) 3 pp_one_piece pp_breaks 1 3 0 0.0 1.0 1.0e-12
  in
  ierr_is "dppgq8 id = 0 ierr" ierr 1;
  near_t 1.0e-14 "dppgq8 integral of 1+2x+3x^2" ans 3.0

let ppgq8_spline_derivatives () =
  let ans, ierr, _ =
    dppgq8 (fun _ -> 1.0) 3 pp_one_piece pp_breaks 1 3 1 0.0 1.0 1.0e-12
  in
  ierr_is "dppgq8 id = 1 ierr" ierr 1;
  near_t 1.0e-14 "dppgq8 integral of 2+6x" ans 5.0;
  let ans, ierr, _ =
    dppgq8 (fun _ -> 1.0) 3 pp_one_piece pp_breaks 1 3 2 0.0 1.0 1.0e-12
  in
  ierr_is "dppgq8 id = 2 ierr" ierr 1;
  near_t 1.0e-14 "dppgq8 integral of 6" ans 6.0

let ppgq8_weighted () =
  let ans, ierr, _ =
    dppgq8 (fun x -> x) 3 pp_one_piece pp_breaks 1 3 0 0.0 1.0 1.0e-12
  in
  ierr_is "dppgq8 weighted ierr" ierr 1;
  near_t 1.0e-14 "dppgq8 integral of x*(1+2x+3x^2)" ans (23.0 /. 12.0)

let ppgq8_reversed () =
  let ans, ierr, _ =
    dppgq8 (fun _ -> 1.0) 3 pp_one_piece pp_breaks 1 3 0 1.0 0.0 1.0e-12
  in
  ierr_is "dppgq8 reversed ierr" ierr 1;
  near_t 1.0e-14 "dppgq8 reversed limits" ans (-3.0)

let ppgq8_error_estimate () =
  let ans, _, e =
    dppgq8 exp 3 pp_one_piece pp_breaks 1 3 0 0.0 1.0 (-1.0e-10)
  in
  let exact = (4.0 *. exp 1.0) -. 5.0 in
  check "dppgq8 error estimate bounds the true error"
    (Float.abs (ans -. exact) <= Float.abs e +. 1.0e-13)

let ppgq8_bad_arguments () =
  let rejects what f =
    match f () with
    | _ -> fail what
    | exception Invalid_argument _ -> ()
  in
  rejects "dppgq8 rejects kk below one" (fun () ->
      dppgq8 (fun _ -> 1.0) 3 pp_one_piece pp_breaks 1 0 0 0.0 1.0 1.0e-12);
  rejects "dppgq8 rejects ldc below kk" (fun () ->
      dppgq8 (fun _ -> 1.0) 2 pp_one_piece pp_breaks 1 3 0 0.0 1.0 1.0e-12);
  rejects "dppgq8 rejects lxi below one" (fun () ->
      dppgq8 (fun _ -> 1.0) 3 pp_one_piece pp_breaks 0 3 0 0.0 1.0 1.0e-12);
  rejects "dppgq8 rejects id at kk" (fun () ->
      dppgq8 (fun _ -> 1.0) 3 pp_one_piece pp_breaks 1 3 3 0.0 1.0 1.0e-12);
  rejects "dppgq8 rejects id below zero" (fun () ->
      dppgq8 (fun _ -> 1.0) 3 pp_one_piece pp_breaks 1 3 (-1) 0.0 1.0 1.0e-12)

(* ---- DPFQAD ---- *)

let pp_two_pieces = [| 1.0; 2.0; 3.0; 4.0 |]
let pp_two_breaks = [| 0.0; 1.0; 2.0 |]

let pfqad_across_breakpoints () =
  let quad, ierr =
    dpfqad (fun _ -> 1.0) 2 pp_two_pieces pp_two_breaks 2 2 0 0.0 2.0 1.0e-12
  in
  ierr_is "dpfqad whole range ierr" ierr 1;
  near_t 1.0e-14 "dpfqad over both pieces" quad 7.0

let pfqad_inside_one_piece () =
  let quad, ierr =
    dpfqad (fun _ -> 1.0) 2 pp_two_pieces pp_two_breaks 2 2 0 0.0 0.5 1.0e-12
  in
  ierr_is "dpfqad one piece ierr" ierr 1;
  near_t 1.0e-14 "dpfqad inside the first piece" quad 0.75

let pfqad_straddling () =
  let quad, ierr =
    dpfqad (fun _ -> 1.0) 2 pp_two_pieces pp_two_breaks 2 2 0 0.5 1.5 1.0e-12
  in
  ierr_is "dpfqad straddling ierr" ierr 1;
  near_t 1.0e-14 "dpfqad straddling a breakpoint" quad 3.25

let pfqad_weighted () =
  let quad, ierr =
    dpfqad (fun x -> x) 2 pp_two_pieces pp_two_breaks 2 2 0 0.0 2.0 1.0e-12
  in
  ierr_is "dpfqad weighted ierr" ierr 1;
  near_t 1.0e-13 "dpfqad weighted by x" quad 9.0

let pfqad_reversed () =
  let quad, ierr =
    dpfqad (fun _ -> 1.0) 2 pp_two_pieces pp_two_breaks 2 2 0 2.0 0.0 1.0e-12
  in
  ierr_is "dpfqad reversed ierr" ierr 1;
  near_t 1.0e-14 "dpfqad reversed limits" quad (-7.0)

let pfqad_equal_limits () =
  let quad, ierr =
    dpfqad (fun _ -> 1.0) 2 pp_two_pieces pp_two_breaks 2 2 0 1.0 1.0 1.0e-12
  in
  check "dpfqad x1 = x2 quad is zero" (quad = 0.0);
  ierr_is "dpfqad x1 = x2 ierr" ierr 1

let pfqad_matches_dppgq8_on_one_piece () =
  let quad, _ =
    dpfqad exp 3 pp_one_piece pp_breaks 1 3 0 0.0 1.0 1.0e-12
  in
  let want, _, _ = dppgq8 exp 3 pp_one_piece pp_breaks 1 3 0 0.0 1.0 1.0e-12 in
  near_t 1.0e-14 "dpfqad matches dppgq8 on a single piece" quad want

let pfqad_derivative_is_a_difference () =
  let quad, ierr =
    dpfqad (fun _ -> 1.0) 2 pp_two_pieces pp_two_breaks 2 2 1 0.0 1.0 1.0e-12
  in
  ierr_is "dpfqad derivative ierr" ierr 1;
  near_t 1.0e-14 "dpfqad integral of the first derivative" quad 2.0

let pfqad_bad_arguments () =
  let rejects what f =
    match f () with
    | _ -> fail what
    | exception Invalid_argument _ -> ()
  in
  rejects "dpfqad rejects k below one" (fun () ->
      dpfqad (fun _ -> 1.0) 2 pp_two_pieces pp_two_breaks 2 0 0 0.0 2.0 1.0e-12);
  rejects "dpfqad rejects ldc below k" (fun () ->
      dpfqad (fun _ -> 1.0) 1 pp_two_pieces pp_two_breaks 2 2 0 0.0 2.0 1.0e-12);
  rejects "dpfqad rejects id at k" (fun () ->
      dpfqad (fun _ -> 1.0) 2 pp_two_pieces pp_two_breaks 2 2 2 0.0 2.0 1.0e-12);
  rejects "dpfqad rejects lxi below one" (fun () ->
      dpfqad (fun _ -> 1.0) 2 pp_two_pieces pp_two_breaks 0 2 0 0.0 2.0 1.0e-12);
  rejects "dpfqad rejects tol above 0.1" (fun () ->
      dpfqad (fun _ -> 1.0) 2 pp_two_pieces pp_two_breaks 2 2 0 0.0 2.0 0.5);
  rejects "dpfqad rejects tol below the unit roundoff" (fun () ->
      dpfqad (fun _ -> 1.0) 2 pp_two_pieces pp_two_breaks 2 2 0 0.0 2.0 1.0e-20)

(* ---- The OCaml surface ---- *)

let estimate what (r : Quad.result) =
  match r.Quad.error with
  | Some e -> e
  | None ->
    fail (what ^ " gave no error estimate");
    0.0

let surface_gauss8 () =
  let r = Quad.gauss8 (fun x -> x *. x) 0.0 1.0 in
  near_t 1.0e-12 "surface gauss8 x^2" r.Quad.value (1.0 /. 3.0);
  check "surface gauss8 converged" r.Quad.converged;
  let e = estimate "surface gauss8 x^2" r in
  check "surface gauss8 error estimate is not negative" (e >= 0.0);
  check "surface gauss8 error estimate bounds the true error"
    (Float.abs (r.Quad.value -. (1.0 /. 3.0)) <= e +. 1.0e-14);
  let r = Quad.gauss8 ~tol:1.0e-12 (fun x -> 1.0 /. x) 1.0 10.0 in
  near_t 1.0e-10 "surface gauss8 1/x" r.Quad.value (log 10.0);
  check "surface gauss8 1/x converged" r.Quad.converged

let surface_gauss8_reversed () =
  let f x = sin x in
  let up = Quad.gauss8 f 0.0 pi and down = Quad.gauss8 f pi 0.0 in
  near_t 1.0e-12 "surface gauss8 sin over 0 pi" up.Quad.value 2.0;
  near_t 1.0e-12 "surface gauss8 reversed negates" down.Quad.value (-2.0)

let surface_gauss8_no_width () =
  let r = Quad.gauss8 (fun x -> exp x) 3.0 3.0 in
  check "surface gauss8 no width is zero" (r.Quad.value = 0.0);
  check "surface gauss8 no width converged" r.Quad.converged

let surface_gauss8_does_not_converge () =
  let r = Quad.gauss8 ~tol:1.0e-14 (fun x -> 1.0 /. sqrt x) 1.0e-12 1.0 in
  check "surface gauss8 reports the tolerance missed"
    (not r.Quad.converged);
  check "surface gauss8 still hands back a value"
    (Float.abs (r.Quad.value -. 2.0) < 1.0e-3);
  check "surface gauss8 owns up to how far off it is"
    (estimate "surface gauss8 unconverged" r > 1.0e-12)

let surface_gauss8_too_narrow () =
  let narrow what a b =
    match Quad.gauss8 (fun x -> x) a b with
    | _ -> fail what
    | exception Quad.Too_narrow (ra, rb) ->
      if ra <> a || rb <> b then fail (what ^ " reported the wrong limits")
  in
  narrow "surface gauss8 too narrow" 1.0 (1.0 +. 1.0e-15);
  narrow "surface gauss8 too narrow negative" (-1.0) (-1.0 -. 1.0e-15)

let surface_newton_cotes7 () =
  let r = Quad.newton_cotes7 ~tol:1.0e-13 (fun x -> 1.0 /. x) 1.0 10.0 in
  near_t 1.0e-11 "surface nc7 1/x" r.Quad.value (log 10.0);
  check "surface nc7 converged" r.Quad.converged;
  check "surface nc7 gives no error estimate" (r.Quad.error = None);
  let r = Quad.newton_cotes7 (fun x -> x *. x *. x) 0.0 2.0 in
  near_t 1.0e-12 "surface nc7 x^3" r.Quad.value 4.0;
  let r = Quad.newton_cotes7 (fun x -> exp x) 2.0 2.0 in
  check "surface nc7 no width is zero" (r.Quad.value = 0.0);
  check "surface nc7 no width converged" r.Quad.converged;
  (match Quad.newton_cotes7 (fun x -> x) 1.0 (1.0 +. 1.0e-15) with
  | _ -> fail "surface nc7 too narrow"
  | exception Quad.Too_narrow (_, _) -> ())

let surface_rules_agree () =
  let f x = exp (-.(x *. x)) in
  let g = Quad.gauss8 f 0.0 2.0 and n = Quad.newton_cotes7 f 0.0 2.0 in
  near_t 1.0e-11 "surface the two rules agree" g.Quad.value n.Quad.value

let surface_table () =
  let x = Array.init 21 (fun i -> float_of_int i /. 20.0) in
  let y = Array.map (fun t -> t *. t *. t) x in
  near_t 1.0e-12 "surface table of x^3" (Quad.integrate_table x y 0.0 1.0) 0.25;
  near_t 1.0e-12 "surface table over part of the range"
    (Quad.integrate_table x y 0.25 0.75)
    ((0.75 ** 4.0 /. 4.0) -. (0.25 ** 4.0 /. 4.0));
  check "surface table over no width"
    (Quad.integrate_table x y 0.5 0.5 = 0.0);
  near_t 1.0e-12 "surface table of two points by trapezoid"
    (Quad.integrate_table [| 0.0; 2.0 |] [| 1.0; 3.0 |] 0.0 2.0)
    4.0

let surface_table_errors () =
  let x = [| 0.0; 1.0; 2.0; 3.0 |] and y = [| 0.0; 1.0; 2.0; 3.0 |] in
  let rejects what f =
    match f () with
    | _ -> fail what
    | exception Invalid_argument _ -> ()
  in
  rejects "surface table rejects one point" (fun () ->
      Quad.integrate_table [| 1.0 |] [| 1.0 |] 0.0 1.0);
  rejects "surface table rejects fewer values than abscissas" (fun () ->
      Quad.integrate_table x [| 0.0; 1.0 |] 0.0 1.0);
  rejects "surface table rejects reversed limits" (fun () ->
      Quad.integrate_table x y 2.0 1.0);
  rejects "surface table rejects abscissas that do not increase" (fun () ->
      Quad.integrate_table [| 0.0; 2.0; 1.0; 3.0 |] y 0.0 3.0);
  rejects "surface table rejects too few abscissas between the limits"
    (fun () -> Quad.integrate_table x y 1.1 1.9)

let () =
  l1_gaus8_x_squared ();
  l1_gaus8_sin ();
  l1_gaus8_exp_neg ();
  l2_gaus8_lorentzian ();
  gaus8_degree_fifteen ();
  gaus8_reversed ();
  gaus8_equal_limits ();
  gaus8_too_nearly_equal ();
  gaus8_error_estimate ();
  gaus8_tighter_is_better ();
  gaus8_insufficient_accuracy ();
  qnc79_degree_seven ();
  qnc79_closed_forms ();
  qnc79_counts_evaluations ();
  qnc79_reversed ();
  qnc79_equal_limits ();
  qnc79_too_nearly_equal ();
  qnc79_additive ();
  qnc79_insufficient_accuracy ();
  davint_parabola_exact ();
  davint_interior_limits ();
  davint_uneven_abscissas ();
  davint_trapezoid ();
  davint_equal_limits ();
  davint_bad_bounds ();
  davint_too_few_between ();
  davint_not_increasing ();
  davint_too_few_values ();
  davint_agrees_with_dgaus8 ();
  ppgq8_spline_value ();
  ppgq8_spline_derivatives ();
  ppgq8_weighted ();
  ppgq8_reversed ();
  ppgq8_error_estimate ();
  ppgq8_bad_arguments ();
  pfqad_across_breakpoints ();
  pfqad_inside_one_piece ();
  pfqad_straddling ();
  pfqad_weighted ();
  pfqad_reversed ();
  pfqad_equal_limits ();
  pfqad_matches_dppgq8_on_one_piece ();
  pfqad_derivative_is_a_difference ();
  pfqad_bad_arguments ();
  surface_gauss8 ();
  surface_gauss8_reversed ();
  surface_gauss8_no_width ();
  surface_gauss8_does_not_converge ();
  surface_gauss8_too_narrow ();
  surface_newton_cotes7 ();
  surface_rules_agree ();
  surface_table ();
  surface_table_errors ();
  if !bad = 0 then print_string "quad: all checks passed\n" else exit 1
