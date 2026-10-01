(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Level 1, level 2 and level 3 checks for the LAPACK auxiliary routines,
   with the expected values taken from the reference Fortran specifications
   and from exact arithmetic. *)

open Reeve.Laux

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
let info_is what (a : int) b = if a <> b then fail what

let mat rows =
  let n = Array.length rows in
  let m = Array.length rows.(0) in
  let a = Array.make (n * m) 0.0 in
  Array.iteri (fun i r -> Array.iteri (fun j v -> a.(i + (j * n)) <- v) r) rows;
  a

let ratio what a b = near what (a /. b) 1.0

let raises what f =
  if not (try f (); false with Invalid_argument _ -> true) then fail what

let l1_lamch_exact () =
  near "l1 dlamch base" (dlamch Base) 2.0;
  near "l1 dlamch digits" (dlamch Digits) 53.0;
  near "l1 dlamch round" (dlamch Round) 1.0;
  near "l1 dlamch min exponent" (dlamch Min_exponent) (-1021.0);
  near "l1 dlamch max exponent" (dlamch Max_exponent) 1024.0;
  check "l1 dlamch eps" (dlamch Eps = Float.ldexp 1.0 (-53));
  check "l1 dlamch precision" (dlamch Precision = Float.ldexp 1.0 (-52));
  check "l1 dlamch precision is epsilon_float" (dlamch Precision = epsilon_float);
  check "l1 dlamch overflow" (dlamch Overflow = max_float);
  check "l1 dlamch underflow" (dlamch Underflow = min_float);
  check "l1 dlamch safe min" (dlamch Safe_min = min_float)

let l1_lamch_relations () =
  check "l1 dlamch eps is half precision"
    (dlamch Eps = dlamch Precision /. dlamch Base);
  check "l1 dlamch one plus eps rounds away" (1.0 +. dlamch Eps = 1.0);
  check "l1 dlamch one plus precision is distinct"
    (1.0 +. dlamch Precision <> 1.0);
  check "l1 dlamch safe min reciprocal is finite"
    (1.0 /. dlamch Safe_min <= dlamch Overflow);
  check "l1 dlamch underflow is base to min exponent minus one"
    (dlamch Underflow = Float.ldexp 1.0 (-1022))

let l1_isnan () =
  check "l1 disnan nan" (disnan nan);
  check "l1 disnan zero" (not (disnan 0.0));
  check "l1 disnan one" (not (disnan 1.0));
  check "l1 disnan infinity" (not (disnan infinity));
  check "l1 disnan neg infinity" (not (disnan neg_infinity));
  check "l1 disnan min float" (not (disnan min_float));
  check "l1 dlaisnan unequal" (dlaisnan 1.0 2.0);
  check "l1 dlaisnan equal" (not (dlaisnan 1.0 1.0));
  check "l1 dlaisnan nan nan" (dlaisnan nan nan);
  check "l1 dlaisnan zeroes" (not (dlaisnan 0.0 (-0.0)))

let l1_lapy2 () =
  near "l1 dlapy2 3 4" (dlapy2 3.0 4.0) 5.0;
  near "l1 dlapy2 4 3" (dlapy2 4.0 3.0) 5.0;
  near "l1 dlapy2 negative" (dlapy2 (-3.0) (-4.0)) 5.0;
  near "l1 dlapy2 zero zero" (dlapy2 0.0 0.0) 0.0;
  near "l1 dlapy2 x zero" (dlapy2 (-7.0) 0.0) 7.0;
  near "l1 dlapy2 zero y" (dlapy2 0.0 (-7.0)) 7.0;
  near "l1 dlapy2 5 12" (dlapy2 5.0 12.0) 13.0;
  near "l1 dlapy2 8 15" (dlapy2 8.0 15.0) 17.0

let l1_lapy3 () =
  near "l1 dlapy3 2 3 6" (dlapy3 2.0 3.0 6.0) 7.0;
  near "l1 dlapy3 1 2 2" (dlapy3 1.0 2.0 2.0) 3.0;
  near "l1 dlapy3 zeroes" (dlapy3 0.0 0.0 0.0) 0.0;
  near "l1 dlapy3 3 4 0" (dlapy3 3.0 4.0 0.0) 5.0;
  near "l1 dlapy3 signs" (dlapy3 (-4.0) 4.0 (-2.0)) 6.0

let l1_lassq_small () =
  let s, q = dlassq 2 [| 3.0; 4.0 |] 0 1 0.0 0.0 in
  near "l1 dlassq 3 4 scale" s 1.0;
  near "l1 dlassq 3 4 sumsq" q 25.0;
  let s, q = dlassq 4 [| 1.0; 1.0; 1.0; 1.0 |] 0 1 0.0 0.0 in
  near "l1 dlassq ones scale" s 1.0;
  near "l1 dlassq ones sumsq" q 4.0;
  let s, q = dlassq 0 [| 5.0 |] 0 1 0.0 0.0 in
  near "l1 dlassq empty scale" s 1.0;
  near "l1 dlassq empty sumsq" q 0.0

let l1_lassq_stride () =
  let s, q = dlassq 2 [| 3.0; 99.0; 4.0; 99.0 |] 0 2 0.0 0.0 in
  near "l1 dlassq stride scale" s 1.0;
  near "l1 dlassq stride sumsq" q 25.0;
  let s, q = dlassq 2 [| 3.0; 4.0 |] 0 (-1) 0.0 0.0 in
  near "l1 dlassq reverse scale" s 1.0;
  near "l1 dlassq reverse sumsq" q 25.0;
  let s, q = dlassq 2 [| 9.0; 3.0; 4.0 |] 1 1 0.0 0.0 in
  near "l1 dlassq offset scale" s 1.0;
  near "l1 dlassq offset sumsq" q 25.0

let l1_laswp_small () =
  let a = mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |] in
  dlaswp 2 a 0 2 0 1 [| 1; 1 |] 0 1;
  same "l1 dlaswp swap two rows" a (mat [| [| 3.0; 4.0 |]; [| 1.0; 2.0 |] |]);
  let a = mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |] in
  dlaswp 2 a 0 2 0 1 [| 0; 1 |] 0 1;
  same "l1 dlaswp identity pivots" a
    (mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |]);
  let a = mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |] in
  dlaswp 2 a 0 2 0 1 [| 1; 1 |] 0 0;
  same "l1 dlaswp zero increment" a
    (mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |])

let l1_laset_full () =
  let a = Array.make 9 1.0 in
  dlaset Full 3 3 7.0 5.0 a 0 3;
  same "l1 dlaset full" a
    [| 5.0; 7.0; 7.0; 7.0; 5.0; 7.0; 7.0; 7.0; 5.0 |]

let l1_lacpy_full () =
  let a = mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |] in
  let b = Array.make 4 0.0 in
  dlacpy Full 2 2 a 0 2 b 0 2;
  same "l1 dlacpy full" b a

let l1_lange_small () =
  let a = mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |] in
  let work = Array.make 2 0.0 in
  near "l1 dlange max abs" (dlange Max_abs 2 2 a 0 2 work) 4.0;
  near "l1 dlange one" (dlange One_norm 2 2 a 0 2 work) 6.0;
  near "l1 dlange inf" (dlange Inf_norm 2 2 a 0 2 work) 7.0;
  near "l1 dlange frobenius" (dlange Frobenius 2 2 a 0 2 work) (sqrt 30.0);
  near "l1 dlange empty rows" (dlange Max_abs 0 2 a 0 2 work) 0.0;
  near "l1 dlange empty cols" (dlange One_norm 2 0 a 0 2 work) 0.0

let l1_ieeeck () =
  info_is "l1 ieeeck infinity" (ieeeck 0 0.0 1.0) 1;
  info_is "l1 ieeeck nan" (ieeeck 1 0.0 1.0) 1

let l1_xerbla () =
  raises "l1 xerbla raises" (fun () -> xerbla "dgetrf" 4);
  (match xerbla "dgetrf" 4 with
  | () -> fail "l1 xerbla returned"
  | exception Invalid_argument m ->
    check "l1 xerbla names the routine"
      (String.length m > 6 && String.equal (String.sub m 0 6) "dgetrf");
    check "l1 xerbla names the position"
      (String.length m > 0 && String.contains m '4'))

let l2_lapy2_scaling () =
  ratio "l2 dlapy2 huge" (dlapy2 3.0e200 4.0e200) 5.0e200;
  ratio "l2 dlapy2 tiny" (dlapy2 3.0e-200 4.0e-200) 5.0e-200;
  ratio "l2 dlapy2 mixed" (dlapy2 3.0e-200 4.0e200) 4.0e200;
  check "l2 dlapy2 infinity x" (dlapy2 infinity 1.0 = infinity);
  check "l2 dlapy2 infinity y" (dlapy2 1.0 neg_infinity = infinity);
  check "l2 dlapy2 nan x" (disnan (dlapy2 nan 1.0));
  check "l2 dlapy2 nan y" (disnan (dlapy2 1.0 nan));
  check "l2 dlapy2 nan both" (disnan (dlapy2 nan nan))

let l2_lapy3_scaling () =
  ratio "l2 dlapy3 huge" (dlapy3 2.0e200 3.0e200 6.0e200) 7.0e200;
  ratio "l2 dlapy3 tiny" (dlapy3 2.0e-200 3.0e-200 6.0e-200) 7.0e-200;
  ratio "l2 dlapy3 mixed" (dlapy3 1.0 2.0e200 2.0e200)
    (2.0e200 *. sqrt 2.0);
  check "l2 dlapy3 infinity" (dlapy3 1.0 1.0 infinity = infinity);
  check "l2 dlapy3 nan" (disnan (dlapy3 nan 0.0 0.0))

let l2_lassq_scaling () =
  let s, q = dlassq 2 [| 3.0e200; 4.0e200 |] 0 1 0.0 0.0 in
  ratio "l2 dlassq huge" (s *. sqrt q) 5.0e200;
  let s, q = dlassq 2 [| 3.0e-200; 4.0e-200 |] 0 1 0.0 0.0 in
  ratio "l2 dlassq tiny" (s *. sqrt q) 5.0e-200;
  let s, q = dlassq 2 [| 3.0e200; 4.0e-200 |] 0 1 0.0 0.0 in
  ratio "l2 dlassq big and small" (s *. sqrt q) 3.0e200;
  let s, q = dlassq 2 [| 1.0; 1.0e-200 |] 0 1 0.0 0.0 in
  near "l2 dlassq mid and small scale" s 1.0;
  near "l2 dlassq mid and small sumsq" q 1.0;
  let s, q = dlassq 2 [| 1.0; 3.0e200 |] 0 1 0.0 0.0 in
  ratio "l2 dlassq mid and big" (s *. sqrt q) 3.0e200

let l2_lassq_accumulate () =
  let s, q = dlassq 1 [| 4.0 |] 0 1 1.0 9.0 in
  near "l2 dlassq accumulate scale" s 1.0;
  near "l2 dlassq accumulate sumsq" q 25.0;
  let s, q = dlassq 0 [| 0.0 |] 0 1 2.0 9.0 in
  near "l2 dlassq carry scale" s 2.0;
  near "l2 dlassq carry sumsq" q 9.0;
  let s, q = dlassq 2 [| 3.0; 4.0 |] 0 1 0.0 5.0 in
  near "l2 dlassq zero scale drops sumsq" s 1.0;
  near "l2 dlassq zero scale sumsq" q 25.0;
  let s, q = dlassq 1 [| 2.0 |] 0 1 3.0 4.0 in
  near "l2 dlassq scaled carry scale" s 1.0;
  near "l2 dlassq scaled carry sumsq" q 40.0;
  let s, q = dlassq 1 [| 1.0 |] 0 1 nan 1.0 in
  check "l2 dlassq nan scale passes through" (disnan s);
  near "l2 dlassq nan scale keeps sumsq" q 1.0;
  let s, q = dlassq 1 [| 1.0 |] 0 1 1.0 nan in
  near "l2 dlassq nan sumsq keeps scale" s 1.0;
  check "l2 dlassq nan sumsq passes through" (disnan q)

let l2_laswp_blocked () =
  let n = 40 in
  let a = Array.make (2 * n) 0.0 in
  for j = 0 to n - 1 do
    a.(0 + (j * 2)) <- float_of_int j;
    a.(1 + (j * 2)) <- float_of_int (100 + j)
  done;
  dlaswp n a 0 2 0 1 [| 1; 1 |] 0 1;
  let want = Array.make (2 * n) 0.0 in
  for j = 0 to n - 1 do
    want.(0 + (j * 2)) <- float_of_int (100 + j);
    want.(1 + (j * 2)) <- float_of_int j
  done;
  same "l2 dlaswp forty columns" a want

let l2_laswp_reverse () =
  let a = mat [| [| 1.0 |]; [| 2.0 |]; [| 3.0 |] |] in
  dlaswp 1 a 0 3 0 2 [| 1; 2; 2 |] 0 1;
  same "l2 dlaswp forward chain" a (mat [| [| 2.0 |]; [| 3.0 |]; [| 1.0 |] |]);
  let a = mat [| [| 2.0 |]; [| 3.0 |]; [| 1.0 |] |] in
  dlaswp 1 a 0 3 0 2 [| 1; 2; 2 |] 0 (-1);
  same "l2 dlaswp reverse undoes forward" a
    (mat [| [| 1.0 |]; [| 2.0 |]; [| 3.0 |] |])

let l2_laswp_offset () =
  let a = Array.make 8 0.0 in
  a.(2) <- 1.0;
  a.(3) <- 3.0;
  a.(4) <- 2.0;
  a.(5) <- 4.0;
  dlaswp 2 a 2 2 0 1 [| 1; 1 |] 0 1;
  same "l2 dlaswp base offset" a
    [| 0.0; 0.0; 3.0; 1.0; 4.0; 2.0; 0.0; 0.0 |]

let l2_laswp_pivot_offset () =
  let a = mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |] in
  dlaswp 2 a 0 2 0 1 [| 9; 9; 9; 1; 1 |] 3 1;
  same "l2 dlaswp pivot base offset" a
    (mat [| [| 3.0; 4.0 |]; [| 1.0; 2.0 |] |])

let l2_laset_triangles () =
  let a = Array.make 15 1.0 in
  dlaset Upper 3 3 7.0 5.0 a 0 5;
  same "l2 dlaset upper with padding" a
    [| 5.0; 1.0; 1.0; 1.0; 1.0; 7.0; 5.0; 1.0; 1.0; 1.0; 7.0; 7.0; 5.0; 1.0;
       1.0 |];
  let a = Array.make 15 1.0 in
  dlaset Lower 3 3 7.0 5.0 a 0 5;
  same "l2 dlaset lower with padding" a
    [| 5.0; 7.0; 7.0; 1.0; 1.0; 1.0; 5.0; 7.0; 1.0; 1.0; 1.0; 1.0; 5.0; 1.0;
       1.0 |];
  let a = Array.make 15 1.0 in
  dlaset Full 3 3 7.0 5.0 a 0 5;
  same "l2 dlaset full with padding" a
    [| 5.0; 7.0; 7.0; 1.0; 1.0; 7.0; 5.0; 7.0; 1.0; 1.0; 7.0; 7.0; 5.0; 1.0;
       1.0 |]

let l2_laset_trapezoid () =
  let a = Array.make 8 1.0 in
  dlaset Upper 2 4 7.0 5.0 a 0 2;
  same "l2 dlaset wide upper trapezoid" a
    [| 5.0; 1.0; 7.0; 5.0; 7.0; 7.0; 7.0; 7.0 |];
  let a = Array.make 8 1.0 in
  dlaset Lower 4 2 7.0 5.0 a 0 4;
  same "l2 dlaset tall lower trapezoid" a
    [| 5.0; 7.0; 7.0; 7.0; 1.0; 5.0; 7.0; 7.0 |];
  let a = Array.make 6 1.0 in
  dlaset Upper 3 2 0.0 1.0 a 0 3;
  same "l2 dlaset identity block" a [| 1.0; 1.0; 1.0; 0.0; 1.0; 1.0 |]

let l2_laset_offset () =
  let a = Array.make 7 1.0 in
  dlaset Full 2 2 7.0 5.0 a 3 2;
  same "l2 dlaset base offset" a
    [| 1.0; 1.0; 1.0; 5.0; 7.0; 7.0; 5.0 |]

let l2_lacpy_triangles () =
  let a = mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |] in
  let b = Array.make 6 0.0 in
  dlacpy Upper 2 2 a 0 2 b 0 3;
  same "l2 dlacpy upper" b [| 1.0; 0.0; 0.0; 2.0; 4.0; 0.0 |];
  let b = Array.make 6 0.0 in
  dlacpy Lower 2 2 a 0 2 b 0 3;
  same "l2 dlacpy lower" b [| 1.0; 3.0; 0.0; 0.0; 4.0; 0.0 |];
  let b = Array.make 6 0.0 in
  dlacpy Full 2 2 a 0 2 b 0 3;
  same "l2 dlacpy full unequal leading dimensions" b
    [| 1.0; 3.0; 0.0; 2.0; 4.0; 0.0 |]

let l2_lacpy_offset () =
  let a = Array.make 6 0.0 in
  a.(2) <- 1.0;
  a.(3) <- 3.0;
  a.(4) <- 2.0;
  a.(5) <- 4.0;
  let b = Array.make 5 0.0 in
  dlacpy Full 2 2 a 2 2 b 1 2;
  same "l2 dlacpy base offsets" b [| 0.0; 1.0; 3.0; 2.0; 4.0 |];
  let a = mat [| [| 1.0; 2.0; 3.0 |]; [| 4.0; 5.0; 6.0 |] |] in
  let b = Array.make 6 0.0 in
  dlacpy Lower 2 3 a 0 2 b 0 2;
  same "l2 dlacpy wide lower" b [| 1.0; 4.0; 0.0; 5.0; 0.0; 0.0 |]

let l2_lange_padded () =
  let a = Array.make 9 0.0 in
  a.(0) <- 1.0;
  a.(1) <- 3.0;
  a.(2) <- 99.0;
  a.(3) <- 2.0;
  a.(4) <- 4.0;
  a.(5) <- 99.0;
  let work = Array.make 2 0.0 in
  near "l2 dlange padded max abs" (dlange Max_abs 2 2 a 0 3 work) 4.0;
  near "l2 dlange padded one" (dlange One_norm 2 2 a 0 3 work) 6.0;
  near "l2 dlange padded inf" (dlange Inf_norm 2 2 a 0 3 work) 7.0;
  near "l2 dlange padded frobenius"
    (dlange Frobenius 2 2 a 0 3 work)
    (sqrt 30.0)

let l2_lange_nan () =
  let a = mat [| [| 1.0; nan |]; [| 3.0; 4.0 |] |] in
  let work = Array.make 2 0.0 in
  check "l2 dlange nan max abs" (disnan (dlange Max_abs 2 2 a 0 2 work));
  check "l2 dlange nan one" (disnan (dlange One_norm 2 2 a 0 2 work));
  check "l2 dlange nan inf" (disnan (dlange Inf_norm 2 2 a 0 2 work))

let l2_lange_scaled () =
  let a = mat [| [| 3.0e200; 0.0 |]; [| 0.0; 4.0e200 |] |] in
  let work = Array.make 2 0.0 in
  ratio "l2 dlange frobenius huge" (dlange Frobenius 2 2 a 0 2 work) 5.0e200;
  let a = mat [| [| 3.0e-200; 0.0 |]; [| 0.0; 4.0e-200 |] |] in
  ratio "l2 dlange frobenius tiny" (dlange Frobenius 2 2 a 0 2 work) 5.0e-200;
  let a = mat [| [| 1.0; 0.0 |]; [| 0.0; 1.0 |] |] in
  near "l2 dlange frobenius identity" (dlange Frobenius 2 2 a 0 2 work)
    (sqrt 2.0);
  near "l2 dlange one identity" (dlange One_norm 2 2 a 0 2 work) 1.0

let l2_lange_offset () =
  let a = Array.make 7 0.0 in
  a.(3) <- 1.0;
  a.(4) <- 3.0;
  a.(5) <- 2.0;
  a.(6) <- 4.0;
  let work = Array.make 2 0.0 in
  near "l2 dlange base offset inf" (dlange Inf_norm 2 2 a 3 2 work) 7.0;
  near "l2 dlange base offset one" (dlange One_norm 2 2 a 3 2 work) 6.0

let l2_ilaenv_fixed () =
  info_is "l2 ilaenv shifts" (ilaenv 4 "DHSEQR" "EN" 0 0 0 0) 6;
  info_is "l2 ilaenv min column" (ilaenv 5 "DGESVD" "" 0 0 0 0) 2;
  info_is "l2 ilaenv svd crossover" (ilaenv 6 "DGESVD" "" 10 20 0 0) 16;
  info_is "l2 ilaenv svd crossover square" (ilaenv 6 "DGESVD" "" 100 100 0 0)
    160;
  info_is "l2 ilaenv processors" (ilaenv 7 "DHSEQR" "" 0 0 0 0) 1;
  info_is "l2 ilaenv multishift" (ilaenv 8 "DHSEQR" "" 0 0 0 0) 50;
  info_is "l2 ilaenv divide and conquer" (ilaenv 9 "DGELSD" "" 0 0 0 0) 25;
  info_is "l2 ilaenv ieee nan" (ilaenv 10 "DGELSD" "" 1 0 0 0) 1;
  info_is "l2 ilaenv ieee infinity" (ilaenv 11 "DGELSD" "" 1 0 0 0) 1;
  info_is "l2 ilaenv below range" (ilaenv 0 "DGETRF" "" 0 0 0 0) (-1);
  info_is "l2 ilaenv above range" (ilaenv 18 "DGETRF" "" 0 0 0 0) (-1)

let l2_ilaenv_names () =
  info_is "l2 ilaenv dgetrf" (ilaenv 1 "DGETRF" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv lower case" (ilaenv 1 "dgetrf" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv mixed case stays mixed"
    (ilaenv 1 "Dgetrf" "" 64 64 (-1) (-1))
    1;
  info_is "l2 ilaenv unknown prefix" (ilaenv 1 "XGETRF" "" 64 64 (-1) (-1)) 1;
  info_is "l2 ilaenv unknown prefix min block"
    (ilaenv 2 "XGEQRF" "" 64 64 (-1) (-1))
    1;
  info_is "l2 ilaenv unknown prefix crossover"
    (ilaenv 3 "XGEQRF" "" 64 64 (-1) (-1))
    1;
  info_is "l2 ilaenv zgetrf" (ilaenv 1 "ZGETRF" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv sgetrf" (ilaenv 1 "SGETRF" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv cgetrf" (ilaenv 1 "CGETRF" "" 64 64 (-1) (-1)) 64

let l2_ilaenv_blocks () =
  info_is "l2 ilaenv dgeqrf" (ilaenv 1 "DGEQRF" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv dgeqrf min" (ilaenv 2 "DGEQRF" "" 64 64 (-1) (-1)) 2;
  info_is "l2 ilaenv dgeqrf crossover" (ilaenv 3 "DGEQRF" "" 64 64 (-1) (-1))
    128;
  info_is "l2 ilaenv dgelqf" (ilaenv 1 "DGELQF" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv dgehrd" (ilaenv 1 "DGEHRD" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv dgehrd crossover" (ilaenv 3 "DGEHRD" "" 64 64 (-1) (-1))
    128;
  info_is "l2 ilaenv dgebrd" (ilaenv 1 "DGEBRD" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv dgetri" (ilaenv 1 "DGETRI" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv dpotrf" (ilaenv 1 "DPOTRF" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv dtrtri" (ilaenv 1 "DTRTRI" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv dtrevc" (ilaenv 1 "DTREVC" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv dlauum" (ilaenv 1 "DLAUUM" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv dlatrs" (ilaenv 1 "DLATRS" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv dlarft crossover" (ilaenv 3 "DLARFT" "" 64 64 (-1) (-1))
    64;
  info_is "l2 ilaenv dlaorhr" (ilaenv 1 "DLAORHR_GETRFNP" "" 64 64 (-1) (-1))
    32;
  info_is "l2 ilaenv dstebz" (ilaenv 1 "DSTEBZ" "" 64 64 (-1) (-1)) 1;
  info_is "l2 ilaenv dgeqp3rk is unreachable"
    (ilaenv 1 "DGEQP3RK" "" 64 64 (-1) (-1))
    1

let l2_ilaenv_symmetric () =
  info_is "l2 ilaenv dsytrf" (ilaenv 1 "DSYTRF" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv dsytrf min" (ilaenv 2 "DSYTRF" "" 64 64 (-1) (-1)) 8;
  info_is "l2 ilaenv dsytrf two stage"
    (ilaenv 1 "DSYTRF_AA_2STAGE" "" 64 64 (-1) (-1))
    192;
  info_is "l2 ilaenv dsytrd" (ilaenv 1 "DSYTRD" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv dsytrd crossover" (ilaenv 3 "DSYTRD" "" 64 64 (-1) (-1))
    32;
  info_is "l2 ilaenv dsygst" (ilaenv 1 "DSYGST" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv zhetrf" (ilaenv 1 "ZHETRF" "" 64 64 (-1) (-1)) 64;
  info_is "l2 ilaenv zhetrf two stage"
    (ilaenv 1 "ZHETRF_AA_2STAGE" "" 64 64 (-1) (-1))
    192;
  info_is "l2 ilaenv zhetrd" (ilaenv 1 "ZHETRD" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv zhetrd crossover" (ilaenv 3 "ZHETRD" "" 64 64 (-1) (-1))
    32;
  info_is "l2 ilaenv dhetrd is not real"
    (ilaenv 1 "DHETRD" "" 64 64 (-1) (-1))
    1

let l2_ilaenv_orthogonal () =
  info_is "l2 ilaenv dorgqr" (ilaenv 1 "DORGQR" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv dormqr" (ilaenv 1 "DORMQR" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv dorgbr" (ilaenv 1 "DORGBR" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv dorgqr min" (ilaenv 2 "DORGQR" "" 64 64 (-1) (-1)) 2;
  info_is "l2 ilaenv dorgqr crossover" (ilaenv 3 "DORGQR" "" 64 64 (-1) (-1))
    128;
  info_is "l2 ilaenv dormqr has no crossover"
    (ilaenv 3 "DORMQR" "" 64 64 (-1) (-1))
    0;
  info_is "l2 ilaenv zungqr" (ilaenv 1 "ZUNGQR" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv zunmqr" (ilaenv 1 "ZUNMQR" "" 64 64 (-1) (-1)) 32;
  info_is "l2 ilaenv dorgsv is not in the set"
    (ilaenv 1 "DORGSV" "" 64 64 (-1) (-1))
    1;
  info_is "l2 ilaenv zorgqr wrong prefix"
    (ilaenv 1 "ZORGQR" "" 64 64 (-1) (-1))
    1

let l2_iparmq_fixed () =
  info_is "l2 iparmq nmin" (iparmq 12 "DLAQR0" "EN" 100 1 100 0) 75;
  info_is "l2 iparmq nibble" (iparmq 14 "DLAQR0" "EN" 100 1 100 0) 14;
  info_is "l2 iparmq cost" (iparmq 17 "DLAQR0" "EN" 100 1 100 0) 10;
  info_is "l2 iparmq below range" (iparmq 11 "DLAQR0" "EN" 100 1 100 0) (-1);
  info_is "l2 iparmq above range" (iparmq 18 "DLAQR0" "EN" 100 1 100 0) (-1);
  info_is "l2 ilaenv forwards to iparmq"
    (iparmq 12 "DLAQR0" "EN" 100 1 100 0)
    (ilaenv 12 "DLAQR0" "EN" 100 1 100 0)

let l3_lamch_blue () =
  let s, q = dlassq 1 [| 1.0 |] 0 1 0.0 0.0 in
  near "l3 dlassq unit scale" s 1.0;
  near "l3 dlassq unit sumsq" q 1.0;
  let s, q = dlassq 1 [| max_float |] 0 1 0.0 0.0 in
  ratio "l3 dlassq overflow threshold" (s *. sqrt q) max_float;
  let s, q = dlassq 1 [| min_float |] 0 1 0.0 0.0 in
  ratio "l3 dlassq underflow threshold" (s *. sqrt q) min_float;
  let s, q = dlassq 2 [| max_float; max_float |] 0 1 0.0 0.0 in
  ratio "l3 dlassq two overflow thresholds" (s *. sqrt (q /. 2.0)) max_float

let l3_lassq_pythagorean () =
  let x = [| 1.0; 2.0; 2.0; 4.0 |] in
  let s, q = dlassq 4 x 0 1 0.0 0.0 in
  near "l3 dlassq pythagorean quadruple scale" s 1.0;
  near "l3 dlassq pythagorean quadruple sumsq" q 25.0;
  let s, q = dlassq 3 [| 2.0; 3.0; 6.0 |] 0 1 0.0 0.0 in
  near "l3 dlassq triple scale" s 1.0;
  near "l3 dlassq triple sumsq" q 49.0;
  near "l3 dlassq agrees with dlapy3" (s *. sqrt q) (dlapy3 2.0 3.0 6.0);
  let s, q = dlassq 2 [| 3.0e200; 4.0e200 |] 0 1 0.0 0.0 in
  ratio "l3 dlassq agrees with dlapy2 scaled" (s *. sqrt q)
    (dlapy2 3.0e200 4.0e200)

let l3_lange_consistency () =
  let a =
    mat
      [| [| 2.0; -1.0; 0.0 |]; [| -1.0; 2.0; -1.0 |]; [| 0.0; -1.0; 2.0 |] |]
  in
  let work = Array.make 3 0.0 in
  near "l3 dlange laplacian max abs" (dlange Max_abs 3 3 a 0 3 work) 2.0;
  near "l3 dlange laplacian one" (dlange One_norm 3 3 a 0 3 work) 4.0;
  near "l3 dlange laplacian inf" (dlange Inf_norm 3 3 a 0 3 work) 4.0;
  near "l3 dlange laplacian frobenius"
    (dlange Frobenius 3 3 a 0 3 work)
    (sqrt 16.0);
  let a = mat [| [| 1.0; 2.0; 3.0 |]; [| 4.0; 5.0; 6.0 |] |] in
  near "l3 dlange wide one" (dlange One_norm 2 3 a 0 2 work) 9.0;
  near "l3 dlange wide inf" (dlange Inf_norm 2 3 a 0 2 work) 15.0;
  near "l3 dlange wide frobenius"
    (dlange Frobenius 2 3 a 0 2 work)
    (sqrt 91.0);
  near "l3 dlange wide max abs" (dlange Max_abs 2 3 a 0 2 work) 6.0

let l3_lacpy_laswp_roundtrip () =
  let a = mat [| [| 1.0; 2.0; 3.0 |]; [| 4.0; 5.0; 6.0 |]; [| 7.0; 8.0; 9.0 |] |] in
  let b = Array.make 9 0.0 in
  dlacpy Full 3 3 a 0 3 b 0 3;
  same "l3 dlacpy copies the whole matrix" b a;
  dlaswp 3 b 0 3 0 2 [| 2; 1; 2 |] 0 1;
  same "l3 dlaswp cycles rows" b
    (mat [| [| 7.0; 8.0; 9.0 |]; [| 4.0; 5.0; 6.0 |]; [| 1.0; 2.0; 3.0 |] |]);
  dlaswp 3 b 0 3 0 2 [| 2; 1; 2 |] 0 (-1);
  same "l3 dlaswp reverse restores" b a

let l3_laset_lacpy_identity () =
  let a = Array.make 16 9.0 in
  dlaset Full 4 4 0.0 1.0 a 0 4;
  let b = Array.make 16 0.0 in
  dlacpy Upper 4 4 a 0 4 b 0 4;
  same "l3 dlacpy upper of identity" b a;
  let b = Array.make 16 0.0 in
  dlacpy Lower 4 4 a 0 4 b 0 4;
  same "l3 dlacpy lower of identity" b a;
  let work = Array.make 4 0.0 in
  near "l3 dlange identity one" (dlange One_norm 4 4 a 0 4 work) 1.0;
  near "l3 dlange identity inf" (dlange Inf_norm 4 4 a 0 4 work) 1.0;
  near "l3 dlange identity frobenius" (dlange Frobenius 4 4 a 0 4 work) 2.0

let l3_iparmq_shifts () =
  info_is "l3 iparmq shifts nh one" (iparmq 15 "DLAQR0" "EN" 1 1 1 0) 2;
  info_is "l3 iparmq shifts nh twenty nine"
    (iparmq 15 "DLAQR0" "EN" 29 1 29 0)
    2;
  info_is "l3 iparmq shifts nh thirty" (iparmq 15 "DLAQR0" "EN" 30 1 30 0) 4;
  info_is "l3 iparmq shifts nh sixty" (iparmq 15 "DLAQR0" "EN" 60 1 60 0) 10;
  info_is "l3 iparmq shifts nh one fifty"
    (iparmq 15 "DLAQR0" "EN" 150 1 150 0)
    20;
  info_is "l3 iparmq shifts nh five ninety"
    (iparmq 15 "DLAQR0" "EN" 590 1 590 0)
    64;
  info_is "l3 iparmq shifts nh three thousand"
    (iparmq 15 "DLAQR0" "EN" 3000 1 3000 0)
    128;
  info_is "l3 iparmq shifts nh six thousand"
    (iparmq 15 "DLAQR0" "EN" 6000 1 6000 0)
    256

let l3_iparmq_window () =
  info_is "l3 iparmq window small" (iparmq 13 "DLAQR0" "EN" 100 1 100 0) 10;
  info_is "l3 iparmq window at swap point"
    (iparmq 13 "DLAQR0" "EN" 500 1 500 0)
    54;
  info_is "l3 iparmq window large" (iparmq 13 "DLAQR0" "EN" 600 1 600 0) 96;
  info_is "l3 iparmq window huge" (iparmq 13 "DLAQR0" "EN" 6000 1 6000 0) 384

let l3_iparmq_accumulate () =
  info_is "l3 iparmq acc laqr small" (iparmq 16 "DLAQR0" "EN" 100 1 100 0) 0;
  info_is "l3 iparmq acc laqr large" (iparmq 16 "DLAQR0" "EN" 600 1 600 0) 2;
  info_is "l3 iparmq acc hseqr large" (iparmq 16 "DHSEQR" "EN" 600 1 600 0) 2;
  info_is "l3 iparmq acc gghrd small" (iparmq 16 "DGGHRD" "EN" 5 1 5 0) 1;
  info_is "l3 iparmq acc gghrd large" (iparmq 16 "DGGHRD" "EN" 20 1 20 0) 2;
  info_is "l3 iparmq acc gghd3 large" (iparmq 16 "DGGHD3" "EN" 20 1 20 0) 2;
  info_is "l3 iparmq acc laexc small" (iparmq 16 "DLAEXC" "EN" 5 1 5 0) 0;
  info_is "l3 iparmq acc laexc large" (iparmq 16 "DLAEXC" "EN" 20 1 20 0) 2;
  info_is "l3 iparmq acc unknown" (iparmq 16 "DGETRF" "EN" 600 1 600 0) 0;
  info_is "l3 iparmq acc lower case" (iparmq 16 "dgghrd" "EN" 20 1 20 0) 2

let l3_ilaenv_dimensions () =
  info_is "l3 ilaenv dgbtrf narrow" (ilaenv 1 "DGBTRF" "" 0 0 0 64) 1;
  info_is "l3 ilaenv dgbtrf wide" (ilaenv 1 "DGBTRF" "" 0 0 0 65) 32;
  info_is "l3 ilaenv dpbtrf narrow" (ilaenv 1 "DPBTRF" "" 0 64 0 0) 1;
  info_is "l3 ilaenv dpbtrf wide" (ilaenv 1 "DPBTRF" "" 0 65 0 0) 32;
  info_is "l3 ilaenv dtrsyl" (ilaenv 1 "DTRSYL" "" 1000 1000 0 0) 160;
  info_is "l3 ilaenv dtrsyl floor" (ilaenv 1 "DTRSYL" "" 10 10 0 0) 48;
  info_is "l3 ilaenv dtrsyl ceiling" (ilaenv 1 "DTRSYL" "" 100000 100000 0 0)
    240;
  info_is "l3 ilaenv ztrsyl" (ilaenv 1 "ZTRSYL" "" 1000 1000 0 0) 80;
  info_is "l3 ilaenv ztrsyl floor" (ilaenv 1 "ZTRSYL" "" 10 10 0 0) 24;
  info_is "l3 ilaenv ztrsyl ceiling" (ilaenv 1 "ZTRSYL" "" 100000 100000 0 0)
    80

let l3_ilaenv_tall_skinny () =
  info_is "l3 ilaenv dgeqr small" (ilaenv 1 "DGEQR" "" 100 100 1 0) 100;
  info_is "l3 ilaenv dgeqr tall" (ilaenv 1 "DGEQR" "" 8192 10000 1 0) 8192;
  info_is "l3 ilaenv dgeqr huge" (ilaenv 1 "DGEQR" "" 10000 10000 1 0) 3;
  info_is "l3 ilaenv dgeqr other mode" (ilaenv 1 "DGEQR" "" 100 100 2 0) 1;
  info_is "l3 ilaenv dgelq small" (ilaenv 1 "DGELQ" "" 100 100 2 0) 100;
  info_is "l3 ilaenv dgelq huge" (ilaenv 1 "DGELQ" "" 10000 10000 2 0) 3;
  info_is "l3 ilaenv dgelq other mode" (ilaenv 1 "DGELQ" "" 100 100 1 0) 1;
  info_is "l3 ilaenv dgghd3" (ilaenv 1 "DGGHD3" "" 64 64 0 0) 32;
  info_is "l3 ilaenv dgghd3 min" (ilaenv 2 "DGGHD3" "" 64 64 0 0) 2;
  info_is "l3 ilaenv dgghd3 crossover" (ilaenv 3 "DGGHD3" "" 64 64 0 0) 128

let () =
  l1_lamch_exact ();
  l1_lamch_relations ();
  l1_isnan ();
  l1_lapy2 ();
  l1_lapy3 ();
  l1_lassq_small ();
  l1_lassq_stride ();
  l1_laswp_small ();
  l1_laset_full ();
  l1_lacpy_full ();
  l1_lange_small ();
  l1_ieeeck ();
  l1_xerbla ();
  l2_lapy2_scaling ();
  l2_lapy3_scaling ();
  l2_lassq_scaling ();
  l2_lassq_accumulate ();
  l2_laswp_blocked ();
  l2_laswp_reverse ();
  l2_laswp_offset ();
  l2_laswp_pivot_offset ();
  l2_laset_triangles ();
  l2_laset_trapezoid ();
  l2_laset_offset ();
  l2_lacpy_triangles ();
  l2_lacpy_offset ();
  l2_lange_padded ();
  l2_lange_nan ();
  l2_lange_scaled ();
  l2_lange_offset ();
  l2_ilaenv_fixed ();
  l2_ilaenv_names ();
  l2_ilaenv_blocks ();
  l2_ilaenv_symmetric ();
  l2_ilaenv_orthogonal ();
  l2_iparmq_fixed ();
  l3_lamch_blue ();
  l3_lassq_pythagorean ();
  l3_lange_consistency ();
  l3_lacpy_laswp_roundtrip ();
  l3_laset_lacpy_identity ();
  l3_iparmq_shifts ();
  l3_iparmq_window ();
  l3_iparmq_accumulate ();
  l3_ilaenv_dimensions ();
  l3_ilaenv_tall_skinny ();
  if !bad = 0 then print_string "laux: all checks passed\n" else exit 1
