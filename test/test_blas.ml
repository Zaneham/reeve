(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Level 1 and level 3 checks for BLAS level 1. *)

open Reeve.Blas
open Reeve.Blas.Raw

let tol = 1.0e-12
let bad = ref 0

let fail what =
  Printf.printf "FAIL %s\n" what;
  bad := !bad + 1

let check what c = if not c then fail what
let near what a b = if not (Float.abs (a -. b) < tol) then fail what
let exact what (a : float) b = if not (a = b) then fail what

let same what (a : float array) b =
  if Array.length a <> Array.length b then fail what
  else Array.iteri (fun i x -> near (Printf.sprintf "%s[%d]" what i) x b.(i)) a

(* ---- Level 3, IBM 360 baseline ---- *)

let l3_daxpy () =
  let pi = 3.14159265358979 in
  let x = [| 1.0; 2.0; 3.0; 4.0; 5.0 |] in
  let y = Array.make 5 100.0 in
  daxpy 5 pi x 0 1 y 0 1;
  same "l3 daxpy" y (Array.init 5 (fun i -> 100.0 +. pi *. float_of_int (i + 1)))

let l3_drotg () =
  let (r, _, c, s) = drotg 3.0 4.0 in
  near "l3 drotg 345 r" r 5.0;
  near "l3 drotg 345 c" c 0.6;
  near "l3 drotg 345 s" s 0.8;
  let (r, _, c, s) = drotg 5.0 12.0 in
  near "l3 drotg 51213 r" r 13.0;
  near "l3 drotg 51213 c" c (5.0 /. 13.0);
  near "l3 drotg 51213 s" s (12.0 /. 13.0)

(* ---- Level 1, the strided paths ---- *)

let l1_daxpy_neg () =
  let x = [| 1.0; 2.0; 3.0 |] in
  let y = [| 10.0; 20.0; 30.0 |] in
  daxpy 3 2.0 x 0 (-1) y 0 1;
  same "daxpy incx -1" y [| 16.0; 24.0; 32.0 |]

let l1_daxpy_zero () =
  let x = [| Float.nan; 1.0 |] in
  let y = [| 7.0; 8.0 |] in
  daxpy 2 0.0 x 0 1 y 0 1;
  same "daxpy alpha 0" y [| 7.0; 8.0 |]

let l1_daxpy_offset () =
  let x = [| 9.0; 1.0; 2.0 |] and y = [| 9.0; 9.0; 10.0; 20.0 |] in
  daxpy 2 3.0 x 1 1 y 2 1;
  same "daxpy offsets" y [| 9.0; 9.0; 13.0; 26.0 |];
  same "daxpy leaves x" x [| 9.0; 1.0; 2.0 |]

let l1_dscal_stride () =
  let x = [| 1.0; 2.0; 3.0; 4.0; 5.0; 6.0 |] in
  dscal 3 10.0 x 0 2;
  same "dscal incx 2" x [| 10.0; 2.0; 30.0; 4.0; 50.0; 6.0 |]

let l1_dscal_nonpos () =
  let x = [| 1.0; 2.0; 3.0 |] in
  dscal 3 10.0 x 0 (-1);
  same "dscal incx -1 is a no-op" x [| 1.0; 2.0; 3.0 |];
  dscal 3 10.0 x 0 0;
  same "dscal incx 0 is a no-op" x [| 1.0; 2.0; 3.0 |]

let l1_dscal_unit () =
  let x = [| 1.0; 2.0; 3.0 |] in
  dscal 3 1.0 x 0 1;
  same "dscal alpha 1" x [| 1.0; 2.0; 3.0 |]

let l1_dscal_offset () =
  let x = [| 9.0; 1.0; 9.0; 2.0 |] in
  dscal 2 3.0 x 1 2;
  same "dscal offset stride" x [| 9.0; 3.0; 9.0; 6.0 |]

let l1_dcopy_reverse () =
  let x = [| 1.0; 2.0; 3.0 |] in
  let y = Array.make 3 0.0 in
  dcopy 3 x 0 1 y 0 (-1);
  same "dcopy incy -1" y [| 3.0; 2.0; 1.0 |]

let l1_dcopy_offset () =
  let x = [| 9.0; 1.0; 2.0; 3.0 |] and y = Array.make 5 9.0 in
  dcopy 3 x 1 1 y 2 1;
  same "dcopy offsets" y [| 9.0; 9.0; 1.0; 2.0; 3.0 |]

let l1_dswap () =
  let x = [| 1.0; 2.0 |] and y = [| 9.0; 8.0 |] in
  dswap 2 x 0 1 y 0 1;
  same "dswap x" x [| 9.0; 8.0 |];
  same "dswap y" y [| 1.0; 2.0 |]

let l1_dswap_negative () =
  let x = [| 1.0; 2.0; 3.0 |] and y = [| 4.0; 5.0; 6.0 |] in
  dswap 3 x 0 1 y 0 (-1);
  same "dswap negative x" x [| 6.0; 5.0; 4.0 |];
  same "dswap negative y" y [| 3.0; 2.0; 1.0 |]

let l1_dswap_offset () =
  let x = [| 9.0; 1.0; 9.0; 2.0 |] and y = [| 8.0; 3.0; 8.0; 4.0 |] in
  dswap 2 x 1 2 y 1 2;
  same "dswap offset x" x [| 9.0; 3.0; 9.0; 4.0 |];
  same "dswap offset y" y [| 8.0; 1.0; 8.0; 2.0 |]

let l1_drot () =
  let x = [| 1.0 |] and y = [| 0.0 |] in
  drot 1 x 0 1 y 0 1 0.6 0.8;
  near "drot x" x.(0) 0.6;
  near "drot y" y.(0) (-0.8)

let l1_drot_offset () =
  let x = [| 9.0; 1.0 |] and y = [| 9.0; 9.0; 0.0 |] in
  drot 1 x 1 1 y 2 1 0.6 0.8;
  near "drot offset x" x.(1) 0.6;
  near "drot offset y" y.(2) (-0.8)

let l1_ddot () =
  let x = [| 1.0; 2.0; 3.0 |] and y = [| 4.0; 5.0; 6.0 |] in
  exact "ddot 3" (ddot 3 x 0 1 y 0 1) 32.0;
  exact "ddot 0" (ddot 0 x 0 1 y 0 1) 0.0;
  exact "ddot negative n" (ddot (-1) x 0 1 y 0 1) 0.0

let l1_ddot_unrolled () =
  let x = Array.init 13 (fun i -> float_of_int (i + 1)) in
  let y = Array.make 13 2.0 in
  for n = 0 to 13 do
    let want = float_of_int (n * (n + 1)) in
    exact (Printf.sprintf "ddot unrolled %d" n) (ddot n x 0 1 y 0 1) want
  done

let l1_ddot_strided () =
  let x = [| 1.0; 9.0; 2.0; 9.0; 3.0 |] and y = [| 4.0; 5.0; 6.0 |] in
  exact "ddot incx 2" (ddot 3 x 0 2 y 0 1) 32.0;
  exact "ddot incy -1" (ddot 3 x 0 2 y 0 (-1)) 28.0

let l1_ddot_offset () =
  let x = [| 9.0; 1.0; 2.0; 3.0 |] and y = [| 9.0; 9.0; 4.0; 5.0; 6.0 |] in
  exact "ddot offsets" (ddot 3 x 1 1 y 2 1) 32.0

let l1_dasum () =
  let x = [| 1.0; -2.0; 3.0 |] in
  exact "dasum 3" (dasum 3 x 0 1) 6.0;
  exact "dasum 0" (dasum 0 x 0 1) 0.0;
  exact "dasum incx -1" (dasum 3 x 0 (-1)) 0.0;
  exact "dasum incx 0" (dasum 3 x 0 0) 0.0

let l1_dasum_unrolled () =
  let x = Array.init 15 (fun i -> float_of_int (-(i + 1))) in
  for n = 0 to 15 do
    let want = float_of_int (n * (n + 1) / 2) in
    exact (Printf.sprintf "dasum unrolled %d" n) (dasum n x 0 1) want
  done

let l1_dasum_strided_offset () =
  let x = [| 9.0; -1.0; 9.0; 2.0; 9.0; -3.0 |] in
  exact "dasum offset stride" (dasum 3 x 1 2) 6.0

let l1_dnrm2 () =
  let x = [| 3.0; 4.0 |] in
  exact "dnrm2 345" (dnrm2 2 x 0 1) 5.0;
  exact "dnrm2 0" (dnrm2 0 x 0 1) 0.0;
  exact "dnrm2 1" (dnrm2 1 [| -7.0 |] 0 1) 7.0;
  exact "dnrm2 zeros" (dnrm2 3 (Array.make 3 0.0) 0 1) 0.0

let l1_dnrm2_scaled () =
  let small = Float.ldexp 1.0 (-600) and big = Float.ldexp 1.0 600 in
  exact "dnrm2 underflowing squares"
    (dnrm2 2 [| 3.0 *. small; 4.0 *. small |] 0 1)
    (5.0 *. small);
  exact "dnrm2 overflowing squares"
    (dnrm2 2 [| 3.0 *. big; 4.0 *. big |] 0 1)
    (5.0 *. big);
  exact "dnrm2 mixed magnitudes" (dnrm2 2 [| big; small |] 0 1) big

let l1_dnrm2_strided_offset () =
  let x = [| 9.0; 3.0; 9.0; 4.0 |] in
  exact "dnrm2 offset stride" (dnrm2 2 x 1 2) 5.0;
  exact "dnrm2 incx -1" (dnrm2 2 [| 3.0; 4.0 |] 0 (-1)) 5.0

let l1_idamax () =
  check "idamax first" (idamax 4 [| 1.0; -5.0; 3.0; 5.0 |] 0 1 = 1);
  check "idamax ties keep first" (idamax 3 [| 2.0; 2.0; 1.0 |] 0 1 = 0);
  check "idamax single" (idamax 1 [| -1.0 |] 0 1 = 0);
  check "idamax empty" (idamax 0 [| 1.0 |] 0 1 = -1);
  check "idamax negative n" (idamax (-1) [| 1.0 |] 0 1 = -1);
  check "idamax nonpositive increment" (idamax 2 [| 1.0; 2.0 |] 0 0 = -1)

let l1_idamax_strided_offset () =
  let x = [| 9.0; 1.0; 9.0; -7.0; 9.0; 3.0 |] in
  check "idamax strided" (idamax 3 x 1 2 = 1)

let l1_empty () =
  let x = [| 1.0 |] and y = [| 2.0 |] in
  daxpy 0 3.0 x 0 1 y 0 1;
  dscal 0 3.0 x 0 1;
  dcopy 0 x 0 1 y 0 1;
  dswap 0 x 0 1 y 0 1;
  drot 0 x 0 1 y 0 1 0.6 0.8;
  same "n 0 leaves x" x [| 1.0 |];
  same "n 0 leaves y" y [| 2.0 |]

let l1_drotg_zero () =
  let (r, z, c, s) = drotg 0.0 0.0 in
  near "drotg 0 r" r 0.0;
  near "drotg 0 z" z 0.0;
  near "drotg 0 c" c 1.0;
  near "drotg 0 s" s 0.0

(* ---- The OCaml surface ---- *)

let raises what f =
  match f () with
  | () -> fail what
  | exception Invalid_argument _ -> ()

let s_axpy () =
  let x = [| 1.0; 2.0; 3.0 |] and y = [| 10.0; 20.0; 30.0 |] in
  axpy 2.0 x y;
  same "axpy" y [| 12.0; 24.0; 36.0 |];
  same "axpy leaves x" x [| 1.0; 2.0; 3.0 |]

let s_axpy_strided () =
  let x = [| 1.0; 99.0; 2.0; 99.0 |] and y = [| 10.0; 20.0 |] in
  axpy ~incx:2 3.0 x y;
  same "axpy strided" y [| 13.0; 26.0 |]

let s_dot () =
  exact "dot" (dot [| 1.0; 2.0; 3.0 |] [| 4.0; 5.0; 6.0 |]) 32.0;
  exact "dot empty" (dot [||] [||]) 0.0;
  exact "dot strided"
    (dot ~incy:2 [| 1.0; 2.0 |] [| 3.0; 9.0; 4.0; 9.0 |]) 11.0

let s_nrm2_asum () =
  exact "nrm2" (nrm2 [| 3.0; 4.0 |]) 5.0;
  exact "nrm2 empty" (nrm2 [||]) 0.0;
  exact "nrm2 strided" (nrm2 ~inc:2 [| 3.0; 9.0; 4.0; 9.0 |]) 5.0;
  exact "asum" (asum [| 1.0; -2.0; 3.0 |]) 6.0;
  exact "asum strided" (asum ~inc:3 [| 1.0; 9.0; 9.0; -2.0 |]) 3.0

let s_scal () =
  let x = [| 1.0; 2.0; 3.0 |] in
  scal 2.0 x;
  same "scal" x [| 2.0; 4.0; 6.0 |];
  let y = [| 1.0; 2.0; 3.0; 4.0 |] in
  scal ~inc:2 (-1.0) y;
  same "scal strided" y [| -1.0; 2.0; -3.0; 4.0 |]

let s_copy_swap () =
  let x = [| 1.0; 2.0 |] and y = [| 8.0; 9.0 |] in
  copy x y;
  same "copy" y [| 1.0; 2.0 |];
  let a = [| 1.0; 2.0 |] and b = [| 3.0; 4.0 |] in
  swap a b;
  same "swap a" a [| 3.0; 4.0 |];
  same "swap b" b [| 1.0; 2.0 |]

let s_rot_rotg () =
  let (r, _, c, s) = rotg 3.0 4.0 in
  near "rotg r" r 5.0;
  let x = [| 3.0 |] and y = [| 4.0 |] in
  rot x y c s;
  near "rot x" x.(0) 5.0;
  near "rot y" y.(0) 0.0

let s_iamax () =
  check "iamax" (iamax [| 1.0; -5.0; 3.0; 5.0 |] = 1);
  check "iamax single" (iamax [| -1.0 |] = 0);
  check "iamax strided" (iamax ~inc:2 [| 1.0; 99.0; -7.0; 99.0 |] = 1)

let s_conformance () =
  raises "axpy length" (fun () -> axpy 1.0 [| 1.0 |] [| 1.0; 2.0 |]);
  raises "dot length" (fun () -> ignore (dot [| 1.0 |] [| 1.0; 2.0 |]));
  raises "copy length" (fun () -> copy [| 1.0; 2.0 |] [| 1.0 |]);
  raises "swap length" (fun () -> swap [| 1.0; 2.0 |] [| 1.0 |]);
  raises "rot length" (fun () -> rot [| 1.0 |] [| 1.0; 2.0 |] 1.0 0.0);
  raises "axpy stride mismatch"
    (fun () -> axpy ~incx:2 1.0 [| 1.0; 2.0; 3.0; 4.0 |] [| 1.0 |]);
  raises "zero increment" (fun () -> ignore (nrm2 ~inc:0 [| 1.0 |]));
  raises "negative increment" (fun () -> ignore (asum ~inc:(-1) [| 1.0 |]));
  raises "scal zero increment" (fun () -> scal ~inc:0 2.0 [| 1.0 |]);
  raises "iamax empty" (fun () -> ignore (iamax [||]))

let s_empty_is_quiet () =
  let x = [||] and y = [||] in
  axpy 1.0 x y;
  scal 2.0 x;
  copy x y;
  swap x y;
  rot x y 0.6 0.8;
  check "empty surface calls are quiet" (Array.length x = 0)

let () =
  l3_daxpy ();
  l3_drotg ();
  l1_daxpy_neg ();
  l1_daxpy_zero ();
  l1_daxpy_offset ();
  l1_dscal_stride ();
  l1_dscal_nonpos ();
  l1_dscal_unit ();
  l1_dscal_offset ();
  l1_dcopy_reverse ();
  l1_dcopy_offset ();
  l1_dswap ();
  l1_dswap_negative ();
  l1_dswap_offset ();
  l1_drot ();
  l1_drot_offset ();
  l1_ddot ();
  l1_ddot_unrolled ();
  l1_ddot_strided ();
  l1_ddot_offset ();
  l1_dasum ();
  l1_dasum_unrolled ();
  l1_dasum_strided_offset ();
  l1_dnrm2 ();
  l1_dnrm2_scaled ();
  l1_dnrm2_strided_offset ();
  l1_idamax ();
  l1_idamax_strided_offset ();
  l1_empty ();
  l1_drotg_zero ();
  s_axpy ();
  s_axpy_strided ();
  s_dot ();
  s_nrm2_asum ();
  s_scal ();
  s_copy_swap ();
  s_rot_rotg ();
  s_iamax ();
  s_conformance ();
  s_empty_is_quiet ();
  if !bad = 0 then print_string "blas: all checks passed\n"
  else exit 1
