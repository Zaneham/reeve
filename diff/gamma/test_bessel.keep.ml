(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Level 1, level 2 and level 3 checks for the Bessel routines, with the
   expected values taken from the slatec-modern Fortran test suites. *)

open Slatec.Bessel

let bad = ref 0

let fail what =
  Printf.printf "FAIL %s\n" what;
  bad := !bad + 1

let near what tol a b = if not (Float.abs (a -. b) < tol) then fail what

let rel_err a b =
  if Float.abs b > 0.0 then Float.abs (a -. b) /. Float.abs b else Float.abs (a -. b)

let close what tol a b = if not (rel_err a b < tol) then fail what

let drift what tol a b =
  if not (rel_err a b < tol || Float.abs (a -. b) < tol) then fail what

let raises what f =
  match f () with
  | _ -> fail what
  | exception Invalid_argument _ -> ()

(* ---- Level 1, regression against the stored SLATEC sequence values ---- *)

let l1_tol = 1.0e-4

let l1_dbesi () =
  let y = Array.make 10 0.0 in
  let nz = dbesi 1.0 0.0 1 5 y in
  if nz <> 0 then fail "l1 dbesi nz";
  let r =
    [| 1.2660658777520082; 0.5651591039924851; 0.1357476697670383;
       0.0221684249243320; 0.0027371202210468 |]
  in
  Array.iteri (fun i v -> drift (Printf.sprintf "l1 dbesi[%d]" i) l1_tol y.(i) v) r

let l1_dbesj () =
  let y = Array.make 10 0.0 in
  let nz = dbesj 1.0 0.0 5 y in
  if nz <> 0 then fail "l1 dbesj nz";
  let r =
    [| 0.7651976865579666; 0.4400505857449335; 0.1149034849319005;
       0.0195633539826684; 0.0024766389641099 |]
  in
  Array.iteri (fun i v -> drift (Printf.sprintf "l1 dbesj[%d]" i) l1_tol y.(i) v) r

let l1_dbesk () =
  let y = Array.make 10 0.0 in
  let nz = dbesk 1.0 0.0 1 5 y in
  if nz <> 0 then fail "l1 dbesk nz";
  let r =
    [| 0.4210244382407084; 0.6019072301972346; 1.6248388986351775;
       7.1012628247379448; 44.234148814293529 |]
  in
  Array.iteri (fun i v -> drift (Printf.sprintf "l1 dbesk[%d]" i) l1_tol y.(i) v) r

(* ---- Level 2, A&S tables 9.1 and 9.8 ---- *)

let l2_tol = 1.0e-2

let l2_ibess_table98 () =
  let data =
    [| [| 0.5; 1.0634834; 0.2578954; 0.0319850 |];
       [| 1.0; 1.2660659; 0.5651591; 0.1357477 |];
       [| 1.5; 1.6467233; 0.9816664; 0.3378346 |];
       [| 2.0; 2.2795853; 1.5906368; 0.6889484 |];
       [| 2.5; 3.2898391; 2.5167163; 1.2764661 |];
       [| 3.0; 4.8807926; 3.9533700; 2.2452125 |];
       [| 4.0; 11.3019220; 9.7594652; 6.4221894 |];
       [| 5.0; 27.2398718; 24.3356421; 17.5056150 |] |]
  in
  let y = Array.make 3 0.0 in
  let worst = ref 0.0 in
  Array.iter
    (fun row ->
      ignore (dbesi row.(0) 0.0 1 3 y);
      for k = 0 to 2 do
        worst := Float.max !worst (rel_err y.(k) row.(k + 1))
      done)
    data;
  if not (!worst < l2_tol) then fail "l2 I bessel vs A&S 9.8"

let l2_jbess_table91 () =
  let data =
    [| [| 0.5; 0.9384698; 0.2422685; 0.0306040; 0.0025603; 0.0001601 |];
       [| 1.0; 0.7651977; 0.4400506; 0.1149035; 0.0195634; 0.0024766 |];
       [| 2.0; 0.2238908; 0.5767248; 0.3528340; 0.1289433; 0.0339957 |];
       [| 5.0; -0.1775968; -0.3275791; 0.0465651; 0.3648312; 0.3912327 |];
       [| 10.0; -0.2459358; 0.0434728; 0.2546303; 0.0583794; -0.2196178 |] |]
  in
  let y = Array.make 5 0.0 in
  let worst = ref 0.0 in
  Array.iter
    (fun row ->
      ignore (dbesj row.(0) 0.0 5 y);
      for k = 0 to 4 do
        worst := Float.max !worst (rel_err y.(k) row.(k + 1))
      done)
    data;
  if not (!worst < l2_tol) then fail "l2 J bessel vs A&S 9.1"

let l2_kbess_table98 () =
  let data =
    [| [| 0.5; 0.9244191; 1.6564411; 7.5501836 |];
       [| 1.0; 0.4210244; 0.6019072; 1.6248389 |];
       [| 1.5; 0.2138056; 0.2773878; 0.5836560 |];
       [| 2.0; 0.1138939; 0.1398659; 0.2537598 |];
       [| 3.0; 0.0347395; 0.0401564; 0.0615104 |];
       [| 5.0; 0.0036911; 0.0040446; 0.0053089 |] |]
  in
  let y = Array.make 3 0.0 in
  let worst = ref 0.0 in
  Array.iter
    (fun row ->
      ignore (dbesk row.(0) 0.0 1 3 y);
      for k = 0 to 2 do
        worst := Float.max !worst (rel_err y.(k) row.(k + 1))
      done)
    data;
  if not (!worst < l2_tol) then fail "l2 K bessel vs A&S 9.8"

let l2_ibess_special () =
  let y = Array.make 2 0.0 in
  ignore (dbesi 0.0 0.0 1 1 y);
  near "l2 I_0(0) is 1" 1.0e-14 y.(0) 1.0;
  ignore (dbesi 0.0 1.0 1 1 y);
  near "l2 I_1(0) is 0" 1.0e-14 y.(0) 0.0;
  ignore (dbesi 100.0 0.0 2 1 y);
  close "l2 I_0 asymptotic" 0.01 y.(0)
    (1.0 /. sqrt (2.0 *. 3.14159265358979 *. 100.0))

let l2_jbess_special () =
  let y = Array.make 2 0.0 in
  ignore (dbesj 0.0 0.0 1 y);
  near "l2 J_0(0) is 1" 1.0e-14 y.(0) 1.0;
  ignore (dbesj 0.0 1.0 1 y);
  near "l2 J_1(0) is 0" 1.0e-14 y.(0) 0.0;
  ignore (dbesj 2.4048255576957728 0.0 1 y);
  near "l2 J_0 first zero" 1.0e-10 y.(0) 0.0

let l2_kbess_asymptotic () =
  let y = Array.make 1 0.0 in
  ignore (dbesk 50.0 0.0 2 1 y);
  close "l2 K_0 asymptotic" 0.01 y.(0) (sqrt (3.14159265358979 /. (2.0 *. 50.0)))

(* ---- Level 2, mathematical identities ---- *)

let l2_j_recurrence () =
  let y = Array.make 1 0.0 in
  let j nu =
    ignore (dbesj 2.0 nu 1 y);
    y.(0)
  in
  let j0 = j 0.0 and j1 = j 1.0 and j2 = j 2.0 in
  near "l2 J recurrence" 1.0e-6 (j0 +. j2) (2.0 /. 2.0 *. j1)

let l2_ik_wronskian () =
  let y = Array.make 1 0.0 in
  let i nu =
    ignore (dbesi 1.0 nu 1 1 y);
    y.(0)
  in
  let k nu =
    ignore (dbesk 1.0 nu 1 1 y);
    y.(0)
  in
  let i0 = i 0.0 and i1 = i 1.0 and k0 = k 0.0 and k1 = k 1.0 in
  near "l2 I-K wronskian" 0.1 ((i0 *. k1) +. (i1 *. k0)) 1.0

(* ---- Level 3, golden values from A&S tables 9.1 and 9.8 ---- *)

let l3_tol = 1.0e-12

let l3_j () =
  let y = Array.make 2 0.0 in
  ignore (dbesj 1.0e-15 0.0 1 y);
  close "l3 J_0(0)" 1.0e-10 y.(0) 1.0;
  ignore (dbesj 1.0 0.0 1 y);
  close "l3 J_0(1)" l3_tol y.(0) 0.7651976865579666;
  ignore (dbesj 2.0 0.0 1 y);
  close "l3 J_0(2)" l3_tol y.(0) 0.2238907791412357;
  ignore (dbesj 5.0 0.0 1 y);
  close "l3 J_0(5)" l3_tol y.(0) (-0.1775967713143383);
  ignore (dbesj 1.0 0.0 2 y);
  close "l3 J_1(1)" l3_tol y.(1) 0.4400505857449335;
  ignore (dbesj 2.0 0.0 2 y);
  close "l3 J_1(2)" l3_tol y.(1) 0.5767248077568734

let l3_i () =
  let y = Array.make 2 0.0 in
  ignore (dbesi 1.0 0.0 1 1 y);
  close "l3 I_0(1)" l3_tol y.(0) 1.2660658777520084;
  ignore (dbesi 2.0 0.0 1 1 y);
  close "l3 I_0(2)" l3_tol y.(0) 2.2795853023360673;
  ignore (dbesi 1.0 0.0 1 2 y);
  close "l3 I_1(1)" l3_tol y.(1) 0.5651591039924851;
  ignore (dbesi 2.0 0.0 1 2 y);
  close "l3 I_1(2)" l3_tol y.(1) 1.5906368546373291

let l3_k () =
  let y = Array.make 2 0.0 in
  ignore (dbesk 1.0 0.0 1 1 y);
  close "l3 K_0(1)" l3_tol y.(0) 0.4210244382407084;
  ignore (dbesk 2.0 0.0 1 1 y);
  close "l3 K_0(2)" l3_tol y.(0) 0.1138938727495334;
  ignore (dbesk 1.0 0.0 1 2 y);
  close "l3 K_1(1)" l3_tol y.(1) 0.6019072301972346;
  ignore (dbesk 2.0 0.0 1 2 y);
  close "l3 K_1(2)" l3_tol y.(1) 0.1398658818165224

(* ---- Level 3, the historical A&S values at non-zero first order ---- *)

let l3_hist_tol = 1.0e-8

let l3_historical () =
  let y = Array.make 1 0.0 in
  ignore (dbesj 1.0 0.0 1 y);
  close "l3h J_0(1)" l3_hist_tol y.(0) 0.7651976865579666;
  ignore (dbesj 1.0 1.0 1 y);
  close "l3h J_1(1)" l3_hist_tol y.(0) 0.4400505857449335;
  ignore (dbesj 2.0 0.0 1 y);
  close "l3h J_0(2)" l3_hist_tol y.(0) 0.2238907791412357;
  ignore (dbesi 1.0 0.0 1 1 y);
  close "l3h I_0(1)" l3_hist_tol y.(0) 1.2660658777520084;
  ignore (dbesi 1.0 1.0 1 1 y);
  close "l3h I_1(1)" l3_hist_tol y.(0) 0.5651591039924850;
  ignore (dbesk 1.0 0.0 1 1 y);
  close "l3h K_0(1)" l3_hist_tol y.(0) 0.4210244382407084

(* ---- Large argument values verified against mpmath in DEVIATIONS.md ---- *)

let large_argument () =
  let y = Array.make 1 0.0 in
  ignore (dbesj 50.0 0.0 1 y);
  close "dev J_0(50)" l3_tol y.(0) 0.0558123276692521;
  ignore (dbesj 100.0 0.0 1 y);
  close "dev J_0(100)" l3_tol y.(0) 0.0199858503042233;
  ignore (dbesj 500.0 0.0 1 y);
  close "dev J_0(500)" l3_tol y.(0) (-0.0341005568807320);
  ignore (dbesj 1000.0 0.0 1 y);
  close "dev J_0(1000)" l3_tol y.(0) 0.0247866861524202;
  ignore (dbesj 20.0 9.0 1 y);
  close "dev J_9(20)" l3_tol y.(0) 0.1251262546479942

(* ---- Fixed order zero and one, and the scaling identities they define ---- *)

let fixed_order () =
  close "dbesi0 1" l3_tol (dbesi0 1.0) 1.2660658777520084;
  close "dbesi0 2" l3_tol (dbesi0 2.0) 2.2795853023360673;
  close "dbesi1 1" l3_tol (dbesi1 1.0) 0.5651591039924851;
  close "dbesi1 2" l3_tol (dbesi1 2.0) 1.5906368546373291;
  close "dbesk0 1" l3_tol (dbesk0 1.0) 0.4210244382407084;
  close "dbesk0 2" l3_tol (dbesk0 2.0) 0.1138938727495334;
  close "dbesk1 1" l3_tol (dbesk1 1.0) 0.6019072301972346;
  close "dbesk1 2" l3_tol (dbesk1 2.0) 0.1398658818165224;
  close "dbesi0e 1" l3_tol (dbesi0e 1.0 *. exp 1.0) 1.2660658777520084;
  close "dbesi1e 2" l3_tol (dbesi1e 2.0 *. exp 2.0) 1.5906368546373291;
  close "dbesk0e 1" l3_tol (dbesk0e 1.0 *. exp (-1.0)) 0.4210244382407084;
  close "dbesk1e 2" l3_tol (dbesk1e 2.0 *. exp (-2.0)) 0.1398658818165224

(* ---- Domains the Fortran refuses ---- *)

let domains () =
  let y = Array.make 2 0.0 in
  raises "dbesj n < 1" (fun () -> dbesj 1.0 0.0 0 y);
  raises "dbesj x < 0" (fun () -> dbesj (-1.0) 0.0 1 y);
  raises "dbesj alpha < 0" (fun () -> dbesj 1.0 (-1.0) 1 y);
  raises "dbesi n < 1" (fun () -> dbesi 1.0 0.0 1 0 y);
  raises "dbesi kode" (fun () -> dbesi 1.0 0.0 3 1 y);
  raises "dbesi x < 0" (fun () -> dbesi (-1.0) 0.0 1 1 y);
  raises "dbesi alpha < 0" (fun () -> dbesi 1.0 (-1.0) 1 1 y);
  raises "dbesk kode" (fun () -> dbesk 1.0 0.0 0 1 y);
  raises "dbesk fnu < 0" (fun () -> dbesk 1.0 (-1.0) 1 1 y);
  raises "dbesk x <= 0" (fun () -> dbesk 0.0 0.0 1 1 y);
  raises "dbesk n < 1" (fun () -> dbesk 1.0 0.0 1 0 y);
  raises "dbesk0 x <= 0" (fun () -> dbesk0 0.0);
  raises "dbesk0e x <= 0" (fun () -> dbesk0e (-1.0));
  raises "dbesk1 x <= 0" (fun () -> dbesk1 0.0);
  raises "dbesk1e x <= 0" (fun () -> dbesk1e (-1.0))

let () =
  l1_dbesi ();
  l1_dbesj ();
  l1_dbesk ();
  l2_ibess_table98 ();
  l2_jbess_table91 ();
  l2_kbess_table98 ();
  l2_ibess_special ();
  l2_jbess_special ();
  l2_kbess_asymptotic ();
  l2_j_recurrence ();
  l2_ik_wronskian ();
  l3_j ();
  l3_i ();
  l3_k ();
  l3_historical ();
  large_argument ();
  fixed_order ();
  domains ();
  if !bad = 0 then print_string "bessel: all checks passed\n" else exit 1
