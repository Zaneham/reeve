(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Checks for the level 2 and level 3 BLAS kernels, with every option
   combination covered and every expected value either exact in binary
   arithmetic or obtained by undoing the operation with its inverse. *)

open Reeve.Blasmat
open Reeve.Blasmat.Raw

module M = Reeve.Mat

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

let matld rows lda =
  let m = Array.length rows.(0) in
  let a = Array.make (lda * m) 0.0 in
  Array.iteri (fun i r -> Array.iteri (fun j v -> a.(i + (j * lda)) <- v) r) rows;
  a

let mat rows = matld rows (Array.length rows)
let zeros n = Array.make n 0.0

let raises_named what name f =
  let m = try f (); "" with Invalid_argument m -> m in
  let k = String.length name in
  check what (String.length m > k && String.sub m 0 k = name)

let seed = ref 1

let rnd () =
  seed := ((!seed * 1103515245) + 12345) mod 2147483648;
  float_of_int ((!seed mod 2001) - 1000) /. 100.0

let tri p =
  let a = Array.make (p * p) 0.0 in
  for j = 0 to p - 1 do
    for i = 0 to p - 1 do
      a.(i + (j * p)) <-
        (if i = j then 4.0 +. float_of_int i
         else 0.5 /. float_of_int (i + j + 2))
    done
  done;
  a

let a22 = [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |]
let b22 = [| [| 5.0; 6.0 |]; [| 7.0; 8.0 |] |]

let gemm_nn () =
  let a = mat a22 and b = mat b22 and c = zeros 4 in
  dgemm No_trans No_trans 2 2 2 1.0 a 0 2 b 0 2 0.0 c 0 2;
  same "gemm n n" c (mat [| [| 19.0; 22.0 |]; [| 43.0; 50.0 |] |])

let gemm_tn () =
  let a = mat a22 and b = mat b22 and c = zeros 4 in
  dgemm Trans No_trans 2 2 2 1.0 a 0 2 b 0 2 0.0 c 0 2;
  same "gemm t n" c (mat [| [| 26.0; 30.0 |]; [| 38.0; 44.0 |] |])

let gemm_nt () =
  let a = mat a22 and b = mat b22 and c = zeros 4 in
  dgemm No_trans Trans 2 2 2 1.0 a 0 2 b 0 2 0.0 c 0 2;
  same "gemm n t" c (mat [| [| 17.0; 23.0 |]; [| 39.0; 53.0 |] |])

let gemm_tt () =
  let a = mat a22 and b = mat b22 and c = zeros 4 in
  dgemm Trans Trans 2 2 2 1.0 a 0 2 b 0 2 0.0 c 0 2;
  same "gemm t t" c (mat [| [| 23.0; 31.0 |]; [| 34.0; 46.0 |] |])

let gemm_conj_is_trans () =
  let a = mat a22 and b = mat b22 and c = zeros 4 in
  dgemm Conj_trans Conj_trans 2 2 2 1.0 a 0 2 b 0 2 0.0 c 0 2;
  same "gemm c c" c (mat [| [| 23.0; 31.0 |]; [| 34.0; 46.0 |] |])

let gemm_alpha_beta () =
  let a = mat a22 and b = mat b22 in
  let c = mat [| [| 1.0; 1.0 |]; [| 1.0; 1.0 |] |] in
  dgemm No_trans No_trans 2 2 2 2.0 a 0 2 b 0 2 (-1.0) c 0 2;
  same "gemm alpha beta" c (mat [| [| 37.0; 43.0 |]; [| 85.0; 99.0 |] |])

let gemm_k_zero () =
  let a = mat a22 and b = mat b22 in
  let c = mat [| [| 4.0; 8.0 |]; [| 16.0; 32.0 |] |] in
  dgemm No_trans No_trans 2 2 0 1.0 a 0 2 b 0 2 0.5 c 0 2;
  same "gemm k zero" c (mat [| [| 2.0; 4.0 |]; [| 8.0; 16.0 |] |])

let gemm_alpha_zero () =
  let a = mat a22 and b = mat b22 in
  let c = mat [| [| 4.0; 8.0 |]; [| 16.0; 32.0 |] |] in
  dgemm No_trans No_trans 2 2 2 0.0 a 0 2 b 0 2 0.0 c 0 2;
  same "gemm alpha zero" c (zeros 4)

let gemm_beta_one_noop () =
  let a = mat a22 and b = mat b22 in
  let c = mat [| [| 4.0; 8.0 |]; [| 16.0; 32.0 |] |] in
  dgemm No_trans No_trans 2 2 2 0.0 a 0 2 b 0 2 1.0 c 0 2;
  same "gemm beta one noop" c (mat [| [| 4.0; 8.0 |]; [| 16.0; 32.0 |] |])

let gemm_offset () =
  let lda = 5 in
  let a = Array.make (lda * 5) 99.0 in
  let b = Array.make (lda * 5) 99.0 in
  let c = Array.make (lda * 5) 99.0 in
  let o = 1 + lda in
  Array.iteri (fun i r -> Array.iteri (fun j v -> a.(o + i + (j * lda)) <- v) r) a22;
  Array.iteri (fun i r -> Array.iteri (fun j v -> b.(o + i + (j * lda)) <- v) r) b22;
  dgemm No_trans No_trans 2 2 2 1.0 a o lda b o lda 0.0 c o lda;
  near "gemm offset 11" c.(o) 19.0;
  near "gemm offset 21" c.(o + 1) 43.0;
  near "gemm offset 12" c.(o + lda) 22.0;
  near "gemm offset 22" c.(o + lda + 1) 50.0;
  let touched = ref 0 in
  Array.iteri
    (fun i v ->
      let r = (i - o) mod lda and q = (i - o) / lda in
      if not (i >= o && r >= 0 && r < 2 && q >= 0 && q < 2) then
        if v <> 99.0 then touched := !touched + 1)
    c;
  check "gemm offset leaves the rest alone" (!touched = 0)

let gemm_rectangular () =
  let a = mat [| [| 1.0; 2.0; 3.0 |]; [| 4.0; 5.0; 6.0 |] |] in
  let b = mat [| [| 1.0 |]; [| 2.0 |]; [| 3.0 |] |] in
  let c = zeros 2 in
  dgemm No_trans No_trans 2 1 3 1.0 a 0 2 b 0 3 0.0 c 0 2;
  same "gemm rectangular" c [| 14.0; 32.0 |]

let gemv_notrans () =
  let a = mat a22 and x = [| 1.0; 1.0 |] and y = zeros 2 in
  dgemv No_trans 2 2 1.0 a 0 2 x 0 1 0.0 y 0 1;
  same "gemv no trans" y [| 3.0; 7.0 |]

let gemv_trans () =
  let a = mat a22 and x = [| 1.0; 1.0 |] and y = zeros 2 in
  dgemv Trans 2 2 1.0 a 0 2 x 0 1 0.0 y 0 1;
  same "gemv trans" y [| 4.0; 6.0 |]

let gemv_conj_is_trans () =
  let a = mat a22 and x = [| 1.0; 1.0 |] and y = zeros 2 in
  dgemv Conj_trans 2 2 1.0 a 0 2 x 0 1 0.0 y 0 1;
  same "gemv conj trans" y [| 4.0; 6.0 |]

let gemv_beta () =
  let a = mat a22 and x = [| 1.0; 0.0 |] and y = [| 10.0; 20.0 |] in
  dgemv No_trans 2 2 2.0 a 0 2 x 0 1 0.5 y 0 1;
  same "gemv beta" y [| 7.0; 16.0 |]

let gemv_negative_increments () =
  let a = mat a22 in
  let x = [| 2.0; 1.0 |] in
  let y = zeros 2 in
  dgemv No_trans 2 2 1.0 a 0 2 x 0 (-1) 0.0 y 0 (-1);
  same "gemv negative increments" y [| 11.0; 5.0 |]

let gemv_strided () =
  let a = mat a22 in
  let x = [| 1.0; 99.0; 1.0 |] in
  let y = [| 0.0; 99.0; 0.0 |] in
  dgemv No_trans 2 2 1.0 a 0 2 x 0 2 0.0 y 0 2;
  same "gemv strided" y [| 3.0; 99.0; 7.0 |]

let gemv_alpha_zero_beta_one () =
  let a = mat a22 and x = [| 1.0; 1.0 |] and y = [| 5.0; 6.0 |] in
  dgemv No_trans 2 2 0.0 a 0 2 x 0 1 1.0 y 0 1;
  same "gemv alpha zero beta one" y [| 5.0; 6.0 |]

let ger_exact () =
  let a = zeros 4 and x = [| 1.0; 2.0 |] and y = [| 3.0; 4.0 |] in
  dger 2 2 1.0 x 0 1 y 0 1 a 0 2;
  same "ger exact" a (mat [| [| 3.0; 4.0 |]; [| 6.0; 8.0 |] |])

let ger_accumulates () =
  let a = mat [| [| 1.0; 1.0 |]; [| 1.0; 1.0 |] |] in
  let x = [| 1.0; 2.0 |] and y = [| 3.0; 4.0 |] in
  dger 2 2 (-1.0) x 0 1 y 0 1 a 0 2;
  same "ger accumulates" a (mat [| [| -2.0; -3.0 |]; [| -5.0; -7.0 |] |])

let ger_negative_increments () =
  let a = zeros 4 and x = [| 2.0; 1.0 |] and y = [| 4.0; 3.0 |] in
  dger 2 2 1.0 x 0 (-1) y 0 (-1) a 0 2;
  same "ger negative increments" a (mat [| [| 3.0; 4.0 |]; [| 6.0; 8.0 |] |])

let ger_alpha_zero () =
  let a = mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |] in
  let x = [| 1.0; 2.0 |] and y = [| 3.0; 4.0 |] in
  dger 2 2 0.0 x 0 1 y 0 1 a 0 2;
  same "ger alpha zero" a (mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |])

let syrk_upper_notrans () =
  let a = mat a22 and c = Array.make 4 99.0 in
  dsyrk Upper No_trans 2 2 1.0 a 0 2 0.0 c 0 2;
  near "syrk u n 11" c.(0) 5.0;
  near "syrk u n 12" c.(2) 11.0;
  near "syrk u n 22" c.(3) 25.0;
  near "syrk u n keeps lower" c.(1) 99.0

let syrk_lower_notrans () =
  let a = mat a22 and c = Array.make 4 99.0 in
  dsyrk Lower No_trans 2 2 1.0 a 0 2 0.0 c 0 2;
  near "syrk l n 11" c.(0) 5.0;
  near "syrk l n 21" c.(1) 11.0;
  near "syrk l n 22" c.(3) 25.0;
  near "syrk l n keeps upper" c.(2) 99.0

let syrk_upper_trans () =
  let a = mat a22 and c = Array.make 4 99.0 in
  dsyrk Upper Trans 2 2 1.0 a 0 2 0.0 c 0 2;
  near "syrk u t 11" c.(0) 10.0;
  near "syrk u t 12" c.(2) 14.0;
  near "syrk u t 22" c.(3) 20.0;
  near "syrk u t keeps lower" c.(1) 99.0

let syrk_lower_trans () =
  let a = mat a22 and c = Array.make 4 99.0 in
  dsyrk Lower Trans 2 2 1.0 a 0 2 0.0 c 0 2;
  near "syrk l t 11" c.(0) 10.0;
  near "syrk l t 21" c.(1) 14.0;
  near "syrk l t 22" c.(3) 20.0;
  near "syrk l t keeps upper" c.(2) 99.0

let syrk_beta () =
  let a = mat a22 in
  let c = mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |] in
  dsyrk Upper No_trans 2 2 (-1.0) a 0 2 2.0 c 0 2;
  near "syrk beta 11" c.(0) (-3.0);
  near "syrk beta 12" c.(2) (-7.0);
  near "syrk beta 22" c.(3) (-17.0);
  near "syrk beta keeps lower" c.(1) 3.0

let syrk_k_zero () =
  let a = mat a22 in
  let c = mat [| [| 4.0; 8.0 |]; [| 16.0; 32.0 |] |] in
  dsyrk Upper No_trans 2 0 1.0 a 0 2 0.5 c 0 2;
  near "syrk k zero 11" c.(0) 2.0;
  near "syrk k zero 12" c.(2) 4.0;
  near "syrk k zero 22" c.(3) 16.0;
  near "syrk k zero keeps lower" c.(1) 16.0

let syrk_rectangular () =
  let a = mat [| [| 1.0; 2.0; 3.0 |]; [| 4.0; 5.0; 6.0 |] |] in
  let c = Array.make 4 99.0 in
  dsyrk Lower No_trans 2 3 1.0 a 0 2 0.0 c 0 2;
  near "syrk rect 11" c.(0) 14.0;
  near "syrk rect 21" c.(1) 32.0;
  near "syrk rect 22" c.(3) 77.0

let trsm_exact_left_lower () =
  let a = mat [| [| 2.0; 0.0 |]; [| 2.0; 4.0 |] |] in
  let b = mat [| [| 4.0 |]; [| 12.0 |] |] in
  dtrsm Left Lower No_trans Non_unit 2 1 1.0 a 0 2 b 0 2;
  same "trsm left lower" b [| 2.0; 2.0 |]

let trsm_exact_left_upper () =
  let a = mat [| [| 2.0; 4.0 |]; [| 0.0; 4.0 |] |] in
  let b = mat [| [| 12.0 |]; [| 8.0 |] |] in
  dtrsm Left Upper No_trans Non_unit 2 1 1.0 a 0 2 b 0 2;
  same "trsm left upper" b [| 2.0; 2.0 |]

let trsm_exact_unit_diag () =
  let a = mat [| [| 9.0; 0.0 |]; [| 2.0; 9.0 |] |] in
  let b = mat [| [| 1.0 |]; [| 4.0 |] |] in
  dtrsm Left Lower No_trans Unit 2 1 1.0 a 0 2 b 0 2;
  same "trsm unit diag" b [| 1.0; 2.0 |]

let trsm_alpha_zero () =
  let a = mat [| [| 2.0; 0.0 |]; [| 2.0; 4.0 |] |] in
  let b = mat [| [| 4.0 |]; [| 12.0 |] |] in
  dtrsm Left Lower No_trans Non_unit 2 1 0.0 a 0 2 b 0 2;
  same "trsm alpha zero" b [| 0.0; 0.0 |]

let trsm_right_exact () =
  let a = mat [| [| 2.0; 0.0 |]; [| 2.0; 4.0 |] |] in
  let b = mat [| [| 8.0; 8.0 |] |] in
  dtrsm Right Lower No_trans Non_unit 1 2 1.0 a 0 2 b 0 1;
  same "trsm right lower" b [| 2.0; 2.0 |]

let trmm_exact_left_lower () =
  let a = mat [| [| 2.0; 0.0 |]; [| 2.0; 4.0 |] |] in
  let b = mat [| [| 2.0 |]; [| 2.0 |] |] in
  dtrmm Left Lower No_trans Non_unit 2 1 1.0 a 0 2 b 0 2;
  same "trmm left lower" b [| 4.0; 12.0 |]

let trmm_exact_left_upper_trans () =
  let a = mat [| [| 2.0; 4.0 |]; [| 0.0; 4.0 |] |] in
  let b = mat [| [| 1.0 |]; [| 1.0 |] |] in
  dtrmm Left Upper Trans Non_unit 2 1 1.0 a 0 2 b 0 2;
  same "trmm left upper trans" b [| 2.0; 8.0 |]

let trmm_alpha_zero () =
  let a = mat [| [| 2.0; 0.0 |]; [| 2.0; 4.0 |] |] in
  let b = mat [| [| 2.0 |]; [| 2.0 |] |] in
  dtrmm Left Lower No_trans Non_unit 2 1 0.0 a 0 2 b 0 2;
  same "trmm alpha zero" b [| 0.0; 0.0 |]

let sides = [| Left; Right |]
let uplos = [| Upper; Lower |]
let transes = [| No_trans; Trans; Conj_trans |]
let diags = [| Unit; Non_unit |]

let name_of si ui ti di m n =
  Printf.sprintf "%s %s %s %s %dx%d"
    (if si = 0 then "left" else "right")
    (if ui = 0 then "upper" else "lower")
    (if ti = 0 then "notrans" else if ti = 1 then "trans" else "conjtrans")
    (if di = 0 then "unit" else "nonunit")
    m n

let trsm_trmm_roundtrip () =
  seed := 7;
  Array.iteri
    (fun si s ->
      Array.iteri
        (fun ui u ->
          Array.iteri
            (fun ti t ->
              Array.iteri
                (fun di d ->
                  List.iter
                    (fun (m, n) ->
                      let p = if si = 0 then m else n in
                      let a = tri (max p 1) in
                      let b = Array.init (m * n) (fun _ -> rnd ()) in
                      let orig = Array.copy b in
                      dtrsm s u t d m n 1.0 a 0 (max p 1) b 0 (max m 1);
                      dtrmm s u t d m n 1.0 a 0 (max p 1) b 0 (max m 1);
                      same_t 1.0e-9
                        ("trsm then trmm " ^ name_of si ui ti di m n)
                        b orig)
                    [ (1, 1); (1, 3); (3, 1); (3, 3); (5, 2); (2, 5); (5, 5) ])
                diags)
            transes)
        uplos)
    sides

let trmm_trsm_roundtrip () =
  seed := 11;
  Array.iteri
    (fun si s ->
      Array.iteri
        (fun ui u ->
          Array.iteri
            (fun ti t ->
              Array.iteri
                (fun di d ->
                  List.iter
                    (fun (m, n) ->
                      let p = if si = 0 then m else n in
                      let a = tri (max p 1) in
                      let b = Array.init (m * n) (fun _ -> rnd ()) in
                      let orig = Array.copy b in
                      dtrmm s u t d m n 2.0 a 0 (max p 1) b 0 (max m 1);
                      dtrsm s u t d m n 2.0 a 0 (max p 1) b 0 (max m 1);
                      same_t 1.0e-9
                        ("trmm then trsm " ^ name_of si ui ti di m n)
                        b (Array.map (fun v -> 4.0 *. v) orig))
                    [ (1, 1); (3, 3); (5, 2); (2, 5) ])
                diags)
            transes)
        uplos)
    sides

let bad_arguments () =
  raises_named "gemm bad m" "dgemm" (fun () ->
      dgemm No_trans No_trans (-1) 1 1 1.0 [| 0.0 |] 0 1 [| 0.0 |] 0 1 0.0
        [| 0.0 |] 0 1);
  raises_named "gemm bad lda" "dgemm" (fun () ->
      dgemm No_trans No_trans 3 1 1 1.0 [| 0.0 |] 0 1 [| 0.0 |] 0 1 0.0
        (zeros 3) 0 3);
  raises_named "gemm bad ldc" "dgemm" (fun () ->
      dgemm No_trans No_trans 3 1 1 1.0 (zeros 3) 0 3 [| 0.0 |] 0 1 0.0
        [| 0.0 |] 0 1);
  raises_named "gemv bad n" "dgemv" (fun () ->
      dgemv No_trans 1 (-1) 1.0 [| 0.0 |] 0 1 [| 0.0 |] 0 1 0.0 [| 0.0 |] 0 1);
  raises_named "gemv zero incx" "dgemv" (fun () ->
      dgemv No_trans 1 1 1.0 [| 0.0 |] 0 1 [| 0.0 |] 0 0 0.0 [| 0.0 |] 0 1);
  raises_named "gemv zero incy" "dgemv" (fun () ->
      dgemv No_trans 1 1 1.0 [| 0.0 |] 0 1 [| 0.0 |] 0 1 0.0 [| 0.0 |] 0 0);
  raises_named "ger bad m" "dger" (fun () ->
      dger (-1) 1 1.0 [| 0.0 |] 0 1 [| 0.0 |] 0 1 [| 0.0 |] 0 1);
  raises_named "ger zero incy" "dger" (fun () ->
      dger 1 1 1.0 [| 0.0 |] 0 1 [| 0.0 |] 0 0 [| 0.0 |] 0 1);
  raises_named "syrk bad k" "dsyrk" (fun () ->
      dsyrk Upper No_trans 1 (-1) 1.0 [| 0.0 |] 0 1 0.0 [| 0.0 |] 0 1);
  raises_named "syrk bad ldc" "dsyrk" (fun () ->
      dsyrk Upper No_trans 3 1 1.0 (zeros 3) 0 3 0.0 [| 0.0 |] 0 1);
  raises_named "trsm bad n" "dtrsm" (fun () ->
      dtrsm Left Lower No_trans Unit 1 (-1) 1.0 [| 0.0 |] 0 1 [| 0.0 |] 0 1);
  raises_named "trsm bad ldb" "dtrsm" (fun () ->
      dtrsm Left Lower No_trans Unit 3 1 1.0 (zeros 9) 0 3 [| 0.0 |] 0 1);
  raises_named "trmm bad m" "dtrmm" (fun () ->
      dtrmm Left Lower No_trans Unit (-1) 1 1.0 [| 0.0 |] 0 1 [| 0.0 |] 0 1);
  raises_named "trmm bad lda" "dtrmm" (fun () ->
      dtrmm Right Lower No_trans Unit 1 3 1.0 (zeros 3) 0 1 (zeros 3) 0 1)

let zero_dimensions_are_quiet () =
  let a = [| 1.0 |] and b = [| 2.0 |] and c = [| 3.0 |] in
  dgemm No_trans No_trans 0 0 0 1.0 a 0 1 b 0 1 1.0 c 0 1;
  dgemv No_trans 0 0 1.0 a 0 1 b 0 1 1.0 c 0 1;
  dger 0 0 1.0 a 0 1 b 0 1 c 0 1;
  dsyrk Upper No_trans 0 0 1.0 a 0 1 1.0 c 0 1;
  dtrsm Left Upper No_trans Unit 0 0 1.0 a 0 1 c 0 1;
  dtrmm Left Upper No_trans Unit 0 0 1.0 a 0 1 c 0 1;
  same "zero dimensions leave a" a [| 1.0 |];
  same "zero dimensions leave b" b [| 2.0 |];
  same "zero dimensions leave c" c [| 3.0 |]

(* ---- The OCaml surface ---- *)

let mt rows =
  M.init (Array.length rows) (Array.length rows.(0)) (fun i j -> rows.(i).(j))

let mnear what a b =
  check (what ^ " shape") (M.rows a = M.rows b && M.cols a = M.cols b);
  for j = 0 to M.cols a - 1 do
    for i = 0 to M.rows a - 1 do
      near (Printf.sprintf "%s (%d,%d)" what i j) (M.get a i j) (M.get b i j)
    done
  done

let s_gemm () =
  let a = mt a22 and b = mt b22 and c = M.create 2 2 in
  gemm a b c;
  mnear "gemm" c (mt [| [| 19.0; 22.0 |]; [| 43.0; 50.0 |] |])

let s_gemm_trans () =
  let a = mt a22 and b = mt b22 and c = M.create 2 2 in
  gemm ~transa:Trans a b c;
  mnear "gemm transa" c (mt [| [| 26.0; 30.0 |]; [| 38.0; 44.0 |] |]);
  let d = M.create 2 2 in
  gemm ~transb:Trans a b d;
  mnear "gemm transb" d (mt [| [| 17.0; 23.0 |]; [| 39.0; 53.0 |] |])

let s_gemm_scalars () =
  let a = mt a22 and b = mt b22 in
  let c = mt [| [| 1.0; 1.0 |]; [| 1.0; 1.0 |] |] in
  gemm ~alpha:2.0 ~beta:3.0 a b c;
  mnear "gemm scalars" c (mt [| [| 41.0; 47.0 |]; [| 89.0; 103.0 |] |])

let s_gemm_rectangular () =
  let a = mt [| [| 1.0; 2.0; 3.0 |]; [| 4.0; 5.0; 6.0 |] |] in
  let b = mt [| [| 7.0; 8.0 |]; [| 9.0; 10.0 |]; [| 11.0; 12.0 |] |] in
  let c = M.create 2 2 in
  gemm a b c;
  mnear "gemm 2x3 3x2" c (mt [| [| 58.0; 64.0 |]; [| 139.0; 154.0 |] |])

let s_gemv () =
  let a = mt a22 in
  let y = zeros 2 in
  gemv a [| 1.0; 1.0 |] y;
  same "gemv" y [| 3.0; 7.0 |];
  let z = [| 1.0; 1.0 |] in
  gemv ~trans:Trans ~alpha:2.0 ~beta:1.0 a [| 1.0; 1.0 |] z;
  same "gemv trans" z [| 9.0; 13.0 |]

let s_ger () =
  let a = M.create 2 2 in
  ger [| 1.0; 2.0 |] [| 3.0; 4.0 |] a;
  mnear "ger" a (mt [| [| 3.0; 4.0 |]; [| 6.0; 8.0 |] |]);
  ger ~alpha:(-1.0) [| 1.0; 2.0 |] [| 3.0; 4.0 |] a;
  mnear "ger undone" a (M.create 2 2)

let s_syrk () =
  let a = mt a22 and c = M.create 2 2 in
  syrk a c;
  near "syrk 00" (M.get c 0 0) 5.0;
  near "syrk 01" (M.get c 0 1) 11.0;
  near "syrk 11" (M.get c 1 1) 25.0;
  check "syrk leaves the lower triangle" (M.get c 1 0 = 0.0);
  let d = M.create 2 2 in
  syrk ~uplo:Lower ~trans:Trans a d;
  near "syrk lower trans 00" (M.get d 0 0) 10.0;
  near "syrk lower trans 10" (M.get d 1 0) 14.0;
  near "syrk lower trans 11" (M.get d 1 1) 20.0;
  check "syrk leaves the upper triangle" (M.get d 0 1 = 0.0)

let s_trsm_trmm_roundtrip () =
  let p = 5 in
  let a = M.of_array p p (tri p) in
  let b0 = M.init p 3 (fun _ _ -> rnd ()) in
  let b = M.copy b0 in
  trsm ~uplo:Upper ~alpha:2.0 a b;
  trmm ~uplo:Upper ~alpha:0.5 a b;
  mnear "trsm then trmm" b b0;
  let a3 = M.of_array 3 3 (tri 3) in
  let c = M.copy b0 in
  trmm ~side:Right ~uplo:Lower ~trans:Trans ~diag:Non_unit a3 c;
  trsm ~side:Right ~uplo:Lower ~trans:Trans ~diag:Non_unit a3 c;
  mnear "trmm then trsm on the right" c b0;
  let d = M.copy b0 in
  trsm ~uplo:Lower ~diag:Unit a d;
  trmm ~uplo:Lower ~diag:Unit a d;
  mnear "trsm then trmm with a unit diagonal" d b0

let s_conformance () =
  let two = M.create 2 2 and three = M.create 3 3 in
  raises_named "gemm inner" "Blasmat.gemm" (fun () ->
      gemm (M.create 2 3) (M.create 2 2) two);
  raises_named "gemm rows" "Blasmat.gemm" (fun () ->
      gemm (M.create 3 2) (M.create 2 2) two);
  raises_named "gemm cols" "Blasmat.gemm" (fun () ->
      gemm (M.create 2 2) (M.create 2 3) two);
  raises_named "gemv x" "Blasmat.gemv" (fun () -> gemv two (zeros 3) (zeros 2));
  raises_named "gemv y" "Blasmat.gemv" (fun () -> gemv two (zeros 2) (zeros 3));
  raises_named "ger x" "Blasmat.ger" (fun () ->
      ger (zeros 3) (zeros 2) two);
  raises_named "ger y" "Blasmat.ger" (fun () ->
      ger (zeros 2) (zeros 3) two);
  raises_named "syrk c not square" "Blasmat.syrk" (fun () ->
      syrk two (M.create 2 3));
  raises_named "syrk order" "Blasmat.syrk" (fun () -> syrk three two);
  raises_named "trsm a not square" "Blasmat.trsm" (fun () ->
      trsm (M.create 2 3) two);
  raises_named "trsm order" "Blasmat.trsm" (fun () -> trsm three two);
  raises_named "trsm right order" "Blasmat.trsm" (fun () ->
      trsm ~side:Right three two);
  raises_named "trmm a not square" "Blasmat.trmm" (fun () ->
      trmm (M.create 3 2) two);
  raises_named "trmm order" "Blasmat.trmm" (fun () -> trmm three two)

let () =
  gemm_nn ();
  gemm_tn ();
  gemm_nt ();
  gemm_tt ();
  gemm_conj_is_trans ();
  gemm_alpha_beta ();
  gemm_k_zero ();
  gemm_alpha_zero ();
  gemm_beta_one_noop ();
  gemm_offset ();
  gemm_rectangular ();
  gemv_notrans ();
  gemv_trans ();
  gemv_conj_is_trans ();
  gemv_beta ();
  gemv_negative_increments ();
  gemv_strided ();
  gemv_alpha_zero_beta_one ();
  ger_exact ();
  ger_accumulates ();
  ger_negative_increments ();
  ger_alpha_zero ();
  syrk_upper_notrans ();
  syrk_lower_notrans ();
  syrk_upper_trans ();
  syrk_lower_trans ();
  syrk_beta ();
  syrk_k_zero ();
  syrk_rectangular ();
  trsm_exact_left_lower ();
  trsm_exact_left_upper ();
  trsm_exact_unit_diag ();
  trsm_alpha_zero ();
  trsm_right_exact ();
  trmm_exact_left_lower ();
  trmm_exact_left_upper_trans ();
  trmm_alpha_zero ();
  trsm_trmm_roundtrip ();
  trmm_trsm_roundtrip ();
  bad_arguments ();
  zero_dimensions_are_quiet ();
  s_gemm ();
  s_gemm_trans ();
  s_gemm_scalars ();
  s_gemm_rectangular ();
  s_gemv ();
  s_ger ();
  s_syrk ();
  s_trsm_trmm_roundtrip ();
  s_conformance ();
  if !bad = 0 then print_string "blasmat: all checks passed\n" else exit 1
