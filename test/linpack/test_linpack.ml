(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Level 1, level 2 and level 3 checks for the LINPACK routines, with the
   matrices, right hand sides and tolerances taken from the slatec-modern
   Fortran suites. *)

open Reeve.Linpack.Raw

module L = Reeve.Linpack
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
let info_is what (a : int) b = if a <> b then fail what

let mat rows =
  let n = Array.length rows in
  let m = Array.length rows.(0) in
  let a = Array.make (n * m) 0.0 in
  Array.iteri (fun i r -> Array.iteri (fun j v -> a.(i + (j * n)) <- v) r) rows;
  a

let band lda n = Array.make (lda * n) 0.0
let set abd lda i j v = abd.(i - 1 + ((j - 1) * lda)) <- v
let get abd lda i j = abd.(i - 1 + ((j - 1) * lda))
let ones n = Array.make n 1.0

let l1_ge_2x2 () =
  let a = mat [| [| 2.0; 1.0 |]; [| 1.0; 3.0 |] |] in
  let b = [| 5.0; 7.0 |] in
  let ipvt = Array.make 2 0 in
  info_is "l1 dgefa 2x2 info" (dgefa a 2 2 ipvt) 0;
  dgesl a 2 2 ipvt b 0;
  same "l1 dgesl 2x2" b [| 1.6; 1.8 |]

let l1_ge_identity () =
  let a = mat [| [| 1.0; 0.0; 0.0 |]; [| 0.0; 1.0; 0.0 |]; [| 0.0; 0.0; 1.0 |] |] in
  let b = [| 1.0; 2.0; 3.0 |] in
  let ipvt = Array.make 3 0 in
  info_is "l1 dgefa identity info" (dgefa a 3 3 ipvt) 0;
  dgesl a 3 3 ipvt b 0;
  same "l1 dgesl identity" b [| 1.0; 2.0; 3.0 |]

let l1_ge_hilbert () =
  let h i j = 1.0 /. float_of_int (i + j + 1) in
  let a = mat (Array.init 3 (fun i -> Array.init 3 (fun j -> h i j))) in
  let b = Array.init 3 (fun i -> h i 0 +. h i 1 +. h i 2) in
  let ipvt = Array.make 3 0 in
  info_is "l1 dgefa hilbert info" (dgefa a 3 3 ipvt) 0;
  dgesl a 3 3 ipvt b 0;
  same_t 1.0e-8 "l1 dgesl hilbert" b (ones 3)

let l1_po_2x2 () =
  let a = mat [| [| 4.0; 2.0 |]; [| 2.0; 5.0 |] |] in
  let b = [| 6.0; 7.0 |] in
  info_is "l1 dpofa 2x2 info" (dpofa a 2 2) 0;
  dposl a 2 2 b;
  same "l1 dposl 2x2" b (ones 2)

let l1_po_identity () =
  let a = mat [| [| 1.0; 0.0; 0.0 |]; [| 0.0; 1.0; 0.0 |]; [| 0.0; 0.0; 1.0 |] |] in
  let b = [| 1.0; 2.0; 3.0 |] in
  info_is "l1 dpofa identity info" (dpofa a 3 3) 0;
  dposl a 3 3 b;
  same "l1 dposl identity" b [| 1.0; 2.0; 3.0 |]

let l1_po_diagonal () =
  let a = mat [| [| 4.0; 0.0; 0.0 |]; [| 0.0; 9.0; 0.0 |]; [| 0.0; 0.0; 16.0 |] |] in
  let b = [| 8.0; 27.0; 64.0 |] in
  info_is "l1 dpofa diagonal info" (dpofa a 3 3) 0;
  dposl a 3 3 b;
  same "l1 dposl diagonal" b [| 2.0; 3.0; 4.0 |]

let l3_ge_doc_example () =
  let a =
    mat [| [| 1.0; 2.0; 3.0 |]; [| 4.0; 5.0; 6.0 |]; [| 7.0; 8.0; 0.0 |] |]
  in
  let b = [| 14.0; 32.0; 23.0 |] in
  let ipvt = Array.make 3 0 in
  info_is "l3 dgefa doc example info" (dgefa a 3 3 ipvt) 0;
  dgesl a 3 3 ipvt b 0;
  same "l3 dgesl doc example" b [| 1.0; 2.0; 3.0 |]

let l3_ge_pascal () =
  let a =
    mat [| [| 1.0; 1.0; 1.0 |]; [| 1.0; 2.0; 3.0 |]; [| 1.0; 3.0; 6.0 |] |]
  in
  let b = [| 3.0; 6.0; 10.0 |] in
  let ipvt = Array.make 3 0 in
  info_is "l3 dgefa pascal info" (dgefa a 3 3 ipvt) 0;
  dgesl a 3 3 ipvt b 0;
  same "l3 dgesl pascal" b (ones 3)

let l3_ge_vandermonde () =
  let a =
    mat [| [| 1.0; 1.0; 1.0 |]; [| 1.0; 2.0; 4.0 |]; [| 1.0; 3.0; 9.0 |] |]
  in
  let b = [| 6.0; 17.0; 34.0 |] in
  let ipvt = Array.make 3 0 in
  info_is "l3 dgefa vandermonde info" (dgefa a 3 3 ipvt) 0;
  dgesl a 3 3 ipvt b 0;
  same "l3 dgesl vandermonde" b [| 1.0; 2.0; 3.0 |]

let l3_po_lehmer () =
  let a =
    mat
      [| [| 1.0; 0.5; 1.0 /. 3.0 |];
         [| 0.5; 1.0; 2.0 /. 3.0 |];
         [| 1.0 /. 3.0; 2.0 /. 3.0; 1.0 |] |]
  in
  let b = [| 11.0 /. 6.0; 13.0 /. 6.0; 2.0 |] in
  info_is "l3 dpofa lehmer info" (dpofa a 3 3) 0;
  dposl a 3 3 b;
  same "l3 dposl lehmer" b (ones 3)

let l3_po_tridiagonal () =
  let a =
    mat
      [| [| 2.0; -1.0; 0.0 |]; [| -1.0; 2.0; -1.0 |]; [| 0.0; -1.0; 2.0 |] |]
  in
  let b = [| 1.0; 0.0; 1.0 |] in
  info_is "l3 dpofa tridiagonal info" (dpofa a 3 3) 0;
  dposl a 3 3 b;
  same "l3 dposl tridiagonal" b (ones 3)

let l3_po_identity_plus_ones () =
  let a =
    mat [| [| 2.0; 1.0; 1.0 |]; [| 1.0; 2.0; 1.0 |]; [| 1.0; 1.0; 2.0 |] |]
  in
  let b = [| 4.0; 4.0; 4.0 |] in
  info_is "l3 dpofa i plus ones info" (dpofa a 3 3) 0;
  dposl a 3 3 b;
  same "l3 dposl i plus ones" b (ones 3)

let l3_ge_forsythe_moler () =
  let a =
    mat [| [| 0.0001; 1.0; 0.0 |]; [| 1.0; 1.0; 0.0 |]; [| 0.0; 0.0; 1.0 |] |]
  in
  let b = [| 1.0001; 2.0; 1.0 |] in
  let ipvt = Array.make 3 0 in
  info_is "l3 dgefa forsythe moler info" (dgefa a 3 3 ipvt) 0;
  dgesl a 3 3 ipvt b 0;
  same_t 1.0e-8 "l3 dgesl forsythe moler" b (ones 3)

let l3_ge_wilkinson () =
  let a =
    mat [| [| 1.0e-20; 1.0; 0.0 |]; [| 1.0; 1.0; 0.0 |]; [| 0.0; 0.0; 1.0 |] |]
  in
  let b = [| 1.0 +. 1.0e-20; 2.0; 1.0 |] in
  let ipvt = Array.make 3 0 in
  info_is "l3 dgefa wilkinson info" (dgefa a 3 3 ipvt) 0;
  dgesl a 3 3 ipvt b 0;
  same "l3 dgesl wilkinson" b (ones 3)

let l3_ge_upper_triangular () =
  let a =
    mat [| [| 1.0; 2.0; 3.0 |]; [| 0.0; 4.0; 5.0 |]; [| 0.0; 0.0; 6.0 |] |]
  in
  let b = [| 14.0; 23.0; 18.0 |] in
  let ipvt = Array.make 3 0 in
  info_is "l3 dgefa upper triangular info" (dgefa a 3 3 ipvt) 0;
  dgesl a 3 3 ipvt b 0;
  same "l3 dgesl upper triangular" b [| 1.0; 2.0; 3.0 |]

let tridiag_band n ml mu d e =
  let lda = (2 * ml) + mu + 1 in
  let abd = band lda n in
  for i = 1 to n do
    set abd lda (ml + mu + 1) i d
  done;
  for i = 1 to n - 1 do
    set abd lda (ml + mu) (i + 1) e;
    set abd lda (ml + mu + 2) i e
  done;
  abd

let penta_band n ml mu =
  let lda = (2 * ml) + mu + 1 in
  let abd = band lda n in
  for j = 1 to n do
    set abd lda (ml + mu + 1) j 6.0
  done;
  for j = 2 to n do
    set abd lda (ml + mu) j (-2.0)
  done;
  for j = 3 to n do
    set abd lda (ml + mu - 1) j (-1.0)
  done;
  for j = 1 to n - 1 do
    set abd lda (ml + mu + 2) j (-2.0)
  done;
  for j = 1 to n - 2 do
    set abd lda (ml + mu + 3) j (-1.0)
  done;
  abd

let l1_gb_tridiagonal () =
  let n = 4 and ml = 1 and mu = 1 in
  let lda = (2 * ml) + mu + 1 in
  let abd = tridiag_band n ml mu 2.0 (-1.0) in
  let b = [| 1.0; 0.0; 0.0; 1.0 |] in
  let ipvt = Array.make n 0 in
  info_is "l1 dgbfa tridiagonal info" (dgbfa abd lda n ml mu ipvt) 0;
  dgbsl abd lda n ml mu ipvt b 0;
  same "l1 dgbsl tridiagonal" b (ones n)

let l1_gb_pentadiagonal () =
  let n = 5 and ml = 2 and mu = 2 in
  let lda = (2 * ml) + mu + 1 in
  let abd = penta_band n ml mu in
  let b = [| 3.0; 1.0; 0.0; 1.0; 3.0 |] in
  let ipvt = Array.make n 0 in
  info_is "l1 dgbfa pentadiagonal info" (dgbfa abd lda n ml mu ipvt) 0;
  dgbsl abd lda n ml mu ipvt b 0;
  same "l1 dgbsl pentadiagonal" b (ones n)

let l1_gb_diagonal () =
  let n = 3 and ml = 0 and mu = 0 in
  let lda = (2 * ml) + mu + 1 in
  let abd = band lda n in
  set abd lda 1 1 2.0;
  set abd lda 1 2 4.0;
  set abd lda 1 3 8.0;
  let b = [| 4.0; 12.0; 32.0 |] in
  let ipvt = Array.make n 0 in
  info_is "l1 dgbfa diagonal info" (dgbfa abd lda n ml mu ipvt) 0;
  dgbsl abd lda n ml mu ipvt b 0;
  same "l1 dgbsl diagonal" b [| 2.0; 3.0; 4.0 |]

let spd_band n m d e =
  let lda = m + 1 in
  let abd = band lda n in
  for j = 1 to n do
    set abd lda (m + 1) j d
  done;
  for j = 2 to n do
    set abd lda m j e
  done;
  abd

let l1_pb_tridiagonal () =
  let n = 4 and m = 1 in
  let lda = m + 1 in
  let abd = spd_band n m 2.0 (-1.0) in
  let b = [| 1.0; 0.0; 0.0; 1.0 |] in
  info_is "l1 dpbfa tridiagonal info" (dpbfa abd lda n m) 0;
  dpbsl abd lda n m b;
  same "l1 dpbsl tridiagonal" b (ones n)

let l1_pb_diagonal () =
  let n = 3 and m = 0 in
  let lda = m + 1 in
  let abd = band lda n in
  set abd lda 1 1 4.0;
  set abd lda 1 2 9.0;
  set abd lda 1 3 16.0;
  let b = [| 8.0; 27.0; 64.0 |] in
  info_is "l1 dpbfa diagonal info" (dpbfa abd lda n m) 0;
  dpbsl abd lda n m b;
  same "l1 dpbsl diagonal" b [| 2.0; 3.0; 4.0 |]

let l1_pb_bandwidth_two () =
  let n = 4 and m = 2 in
  let lda = m + 1 in
  let abd = band lda n in
  for j = 1 to n do
    set abd lda (m + 1) j 10.0
  done;
  for j = 2 to n do
    set abd lda m j (-3.0)
  done;
  for j = 3 to n do
    set abd lda (m - 1) j 1.0
  done;
  let b = [| 8.0; 5.0; 5.0; 8.0 |] in
  info_is "l1 dpbfa bandwidth two info" (dpbfa abd lda n m) 0;
  dpbsl abd lda n m b;
  same "l1 dpbsl bandwidth two" b (ones n)

let l1_geco_identity () =
  let a = mat [| [| 1.0; 0.0; 0.0 |]; [| 0.0; 1.0; 0.0 |]; [| 0.0; 0.0; 1.0 |] |] in
  let ipvt = Array.make 3 0 and z = Array.make 3 0.0 in
  near_t 0.1 "l1 dgeco identity" (dgeco a 3 3 ipvt z) 1.0

let l1_geco_diagonal () =
  let a = mat [| [| 1.0; 0.0; 0.0 |]; [| 0.0; 2.0; 0.0 |]; [| 0.0; 0.0; 4.0 |] |] in
  let ipvt = Array.make 3 0 and z = Array.make 3 0.0 in
  let r = dgeco a 3 3 ipvt z in
  check "l1 dgeco diagonal" (r > 0.1 && r < 0.5)

let l1_geco_hilbert () =
  let a =
    mat
      [| [| 1.0; 0.5; 1.0 /. 3.0 |];
         [| 0.5; 1.0 /. 3.0; 0.25 |];
         [| 1.0 /. 3.0; 0.25; 0.2 |] |]
  in
  let ipvt = Array.make 3 0 and z = Array.make 3 0.0 in
  let r = dgeco a 3 3 ipvt z in
  check "l1 dgeco hilbert" (r < 0.01 && r > 0.0)

let l2_gb_residual () =
  let n = 4 and ml = 1 and mu = 1 in
  let lda = (2 * ml) + mu + 1 in
  let abd = tridiag_band n ml mu 2.0 (-1.0) in
  let orig = Array.copy abd in
  let b = [| 1.0; 0.0; 0.0; 1.0 |] in
  let x = Array.copy b in
  let ipvt = Array.make n 0 in
  info_is "l2 dgbfa residual info" (dgbfa abd lda n ml mu ipvt) 0;
  dgbsl abd lda n ml mu ipvt x 0;
  let s = ref 0.0 in
  for i = 1 to n do
    let ax = ref (get orig lda (ml + mu + 1) i *. x.(i - 1)) in
    if i > 1 then ax := !ax +. (get orig lda (ml + mu + 2) (i - 1) *. x.(i - 2));
    if i < n then ax := !ax +. (get orig lda (ml + mu) (i + 1) *. x.(i));
    s := !s +. ((!ax -. b.(i - 1)) *. (!ax -. b.(i - 1)))
  done;
  check "l2 dgbsl residual" (sqrt !s < tol)

let l2_gb_linearity () =
  let n = 5 and ml = 2 and mu = 2 in
  let lda = (2 * ml) + mu + 1 in
  let abd = penta_band n ml mu in
  let b1 = [| 3.0; 1.0; 0.0; 1.0; 3.0 |] in
  let b2 = [| 6.0; 2.0; 0.0; 2.0; 6.0 |] in
  let ipvt = Array.make n 0 in
  info_is "l2 dgbfa linearity info" (dgbfa abd lda n ml mu ipvt) 0;
  dgbsl abd lda n ml mu ipvt b1 0;
  dgbsl abd lda n ml mu ipvt b2 0;
  same "l2 dgbsl linearity" b2 (Array.map (fun x -> 2.0 *. x) b1)

let l2_gb_inverse_column () =
  let n = 3 and ml = 1 and mu = 1 in
  let lda = (2 * ml) + mu + 1 in
  let abd = tridiag_band n ml mu 4.0 (-1.0) in
  let c = [| 1.0; 0.0; 0.0 |] in
  let ipvt = Array.make n 0 in
  info_is "l2 dgbfa inverse column info" (dgbfa abd lda n ml mu ipvt) 0;
  dgbsl abd lda n ml mu ipvt c 0;
  check "l2 dgbsl inverse column"
    (Array.for_all (fun v -> Float.abs v < 1.0) c && c.(0) > 0.0)

let l2_gb_zero_rhs () =
  let n = 4 and ml = 1 and mu = 1 in
  let lda = (2 * ml) + mu + 1 in
  let abd = tridiag_band n ml mu 2.0 (-1.0) in
  let b = Array.make n 0.0 in
  let ipvt = Array.make n 0 in
  info_is "l2 dgbfa zero rhs info" (dgbfa abd lda n ml mu ipvt) 0;
  dgbsl abd lda n ml mu ipvt b 0;
  same "l2 dgbsl zero rhs" b (Array.make n 0.0)

let l2_pb_rtr_diagonal () =
  let n = 4 and m = 1 in
  let lda = m + 1 in
  let abd = spd_band n m 2.0 (-1.0) in
  let orig = Array.copy abd in
  info_is "l2 dpbfa rtr info" (dpbfa abd lda n m) 0;
  let err = ref 0.0 in
  for i = 1 to n do
    let d = get abd lda (m + 1) i in
    let e = if i > 1 then get abd lda m i else 0.0 in
    let v = Float.abs (get orig lda (m + 1) i -. ((d *. d) +. (e *. e))) in
    err := Float.max !err v
  done;
  check "l2 dpbfa rtr diagonal" (!err < 1.0e-6)

let l2_pb_residual () =
  let n = 4 and m = 1 in
  let lda = m + 1 in
  let abd = spd_band n m 2.0 (-1.0) in
  let orig = Array.copy abd in
  let b = [| 1.0; 0.0; 0.0; 1.0 |] in
  let x = Array.copy b in
  info_is "l2 dpbfa residual info" (dpbfa abd lda n m) 0;
  dpbsl abd lda n m x;
  let s = ref 0.0 in
  for i = 1 to n do
    let ax = ref (get orig lda (m + 1) i *. x.(i - 1)) in
    if i > 1 then ax := !ax +. (get orig lda m i *. x.(i - 2));
    if i < n then ax := !ax +. (get orig lda m (i + 1) *. x.(i));
    s := !s +. ((!ax -. b.(i - 1)) *. (!ax -. b.(i - 1)))
  done;
  check "l2 dpbsl residual" (sqrt !s < tol)

let l2_pb_inverse_symmetry () =
  let n = 3 and m = 1 in
  let lda = m + 1 in
  let abd = spd_band n m 4.0 (-1.0) in
  info_is "l2 dpbfa inverse symmetry info" (dpbfa abd lda n m) 0;
  let c1 = [| 1.0; 0.0; 0.0 |] and c2 = [| 0.0; 1.0; 0.0 |] in
  dpbsl abd lda n m c1;
  dpbsl abd lda n m c2;
  near "l2 dpbsl inverse symmetry" c1.(1) c2.(0)

let l2_pb_not_definite () =
  let n = 2 and m = 1 in
  let lda = m + 1 in
  let abd = band lda n in
  set abd lda (m + 1) 1 1.0;
  set abd lda (m + 1) 2 1.0;
  set abd lda m 2 2.0;
  check "l2 dpbfa detects indefinite" (dpbfa abd lda n m <> 0)

let l2_geco_range () =
  let a = mat [| [| 1.0; 0.0; 0.0 |]; [| 0.0; 2.0; 0.0 |]; [| 0.0; 0.0; 3.0 |] |] in
  let ipvt = Array.make 3 0 and z = Array.make 3 0.0 in
  let r = dgeco a 3 3 ipvt z in
  check "l2 dgeco rcond range" (r > 0.0 && r <= 1.0)

let l2_geco_scaling () =
  let a1 = mat [| [| 4.0; 1.0 |]; [| 1.0; 3.0 |] |] in
  let a2 = Array.map (fun v -> 10.0 *. v) a1 in
  let ipvt = Array.make 2 0 and z = Array.make 2 0.0 in
  let r1 = dgeco a1 2 2 ipvt z in
  let r2 = dgeco a2 2 2 ipvt z in
  near_t 0.01 "l2 dgeco scaling invariance" r1 r2

let l2_geco_near_singular () =
  let e = epsilon_float in
  let a =
    mat
      [| [| 1.0; 4.0; 7.0 |];
         [| 2.0; 5.0; 8.0 |];
         [| 3.0 +. e; 9.0 +. e; 15.0 +. e |] |]
  in
  let ipvt = Array.make 3 0 and z = Array.make 3 0.0 in
  check "l2 dgeco near singular" (dgeco a 3 3 ipvt z < 1.0e-10)

let l2_geco_orthogonal () =
  let c = cos 0.7 and s = sin 0.7 in
  let a = mat [| [| c; -.s |]; [| s; c |] |] in
  let ipvt = Array.make 2 0 and z = Array.make 2 0.0 in
  check "l2 dgeco orthogonal" (dgeco a 2 2 ipvt z > 0.3)

let l2_geco_triangular () =
  let a =
    mat [| [| 1.0; 1.0; 1.0 |]; [| 0.0; 2.0; 1.0 |]; [| 0.0; 0.0; 4.0 |] |]
  in
  let ipvt = Array.make 3 0 and z = Array.make 3 0.0 in
  let r = dgeco a 3 3 ipvt z in
  check "l2 dgeco triangular" (r > 0.1 && r < 0.5)

let l3_gb_laplacian () =
  let n = 4 and ml = 1 and mu = 1 in
  let lda = (2 * ml) + mu + 1 in
  let h = 0.2 in
  let abd = tridiag_band n ml mu (2.0 /. (h *. h)) (-1.0 /. (h *. h)) in
  let b = Array.make n 2.0 in
  let ipvt = Array.make n 0 in
  info_is "l3 dgbfa laplacian info" (dgbfa abd lda n ml mu ipvt) 0;
  dgbsl abd lda n ml mu ipvt b 0;
  let x =
    Array.init n (fun i ->
        let xi = float_of_int (i + 1) *. h in
        xi *. (1.0 -. xi))
  in
  same_t 1.0e-12 "l3 dgbsl laplacian" b x

let l3_gb_convection_diffusion () =
  let n = 4 and ml = 1 and mu = 1 in
  let lda = (2 * ml) + mu + 1 in
  let h = 0.2 and e = 0.1 in
  let sub = (-.e /. (h *. h)) -. (1.0 /. h) in
  let d = (2.0 *. e /. (h *. h)) +. (1.0 /. h) in
  let sup = -.e /. (h *. h) in
  let abd = band lda n in
  for j = 1 to n do
    set abd lda (ml + mu + 1) j d
  done;
  for j = 2 to n do
    set abd lda (ml + mu) j sup
  done;
  for j = 1 to n - 1 do
    set abd lda (ml + mu + 2) j sub
  done;
  let b = Array.make n 0.0 in
  b.(n - 1) <- -.sup;
  let ipvt = Array.make n 0 in
  info_is "l3 dgbfa convection diffusion info" (dgbfa abd lda n ml mu ipvt) 0;
  dgbsl abd lda n ml mu ipvt b 0;
  check "l3 dgbsl convection diffusion"
    (b.(0) > 0.0 && b.(n - 1) > b.(0) && b.(n - 1) < 1.5)

let l3_gb_laplacian_slice () =
  let n = 3 and ml = 1 and mu = 1 in
  let lda = (2 * ml) + mu + 1 in
  let abd = tridiag_band n ml mu 4.0 (-1.0) in
  let b = [| 3.0; 2.0; 3.0 |] in
  let ipvt = Array.make n 0 in
  info_is "l3 dgbfa laplacian slice info" (dgbfa abd lda n ml mu ipvt) 0;
  dgbsl abd lda n ml mu ipvt b 0;
  same "l3 dgbsl laplacian slice" b (ones n)

let l3_pb_fem_bar () =
  let n = 4 and m = 1 in
  let lda = m + 1 in
  let abd = spd_band n m 2.0 (-1.0) in
  let b = ones n in
  info_is "l3 dpbfa fem bar info" (dpbfa abd lda n m) 0;
  dpbsl abd lda n m b;
  check "l3 dpbsl fem bar positive" (b.(0) > 0.0 && b.(3) > 0.0);
  near "l3 dpbsl fem bar outer symmetry" b.(0) b.(3);
  near "l3 dpbsl fem bar inner symmetry" b.(1) b.(2)

let l3_pb_fem_mass () =
  let n = 3 and m = 1 in
  let lda = m + 1 in
  let abd = spd_band n m (4.0 /. 6.0) (1.0 /. 6.0) in
  let b = ones n in
  info_is "l3 dpbfa fem mass info" (dpbfa abd lda n m) 0;
  dpbsl abd lda n m b;
  check "l3 dpbsl fem mass positive" (Array.for_all (fun v -> v > 0.0) b)

let l3_pb_string_mode () =
  let n = 5 and m = 1 in
  let lda = m + 1 in
  let abd = spd_band n m 2.0 (-1.0) in
  let pi = 4.0 *. atan 1.0 in
  let x =
    Array.init n (fun i ->
        sin (pi *. float_of_int (i + 1) /. float_of_int (n + 1)))
  in
  let b = Array.copy x in
  info_is "l3 dpbfa string mode info" (dpbfa abd lda n m) 0;
  dpbsl abd lda n m b;
  near_t 0.01 "l3 dpbsl string mode" (b.(2) /. x.(2)) (b.(0) /. x.(0))

let l3_geco_hilbert () =
  let a =
    mat
      (Array.init 3 (fun i ->
           Array.init 3 (fun j -> 1.0 /. float_of_int (i + j + 1))))
  in
  let ipvt = Array.make 3 0 and z = Array.make 3 0.0 in
  let r = dgeco a 3 3 ipvt z in
  check "l3 dgeco hilbert" (r < 0.01 && r > 0.0001)

let l3_geco_vandermonde () =
  let nodes = [| 1.0; 2.0; 3.0 |] in
  let a =
    mat
      (Array.init 3 (fun i ->
           Array.init 3 (fun j -> Float.pow nodes.(i) (float_of_int j))))
  in
  let ipvt = Array.make 3 0 and z = Array.make 3 0.0 in
  let r = dgeco a 3 3 ipvt z in
  check "l3 dgeco vandermonde" (r > 0.001 && r < 0.5)

let l3_geco_moler () =
  let a =
    mat
      (Array.init 3 (fun i ->
           Array.init 3 (fun j -> float_of_int (min (i + 1) (j + 1)) -. 2.0)))
  in
  let ipvt = Array.make 3 0 and z = Array.make 3 0.0 in
  let r = dgeco a 3 3 ipvt z in
  check "l3 dgeco moler" (r > 0.0 && r < 1.0)

let l3_geco_pascal () =
  let a =
    mat [| [| 1.0; 1.0; 1.0 |]; [| 1.0; 2.0; 3.0 |]; [| 1.0; 3.0; 6.0 |] |]
  in
  let ipvt = Array.make 3 0 and z = Array.make 3 0.0 in
  check "l3 dgeco pascal" (dgeco a 3 3 ipvt z > 0.01)

let raises what f =
  if not (try f (); false with Invalid_argument _ -> true) then fail what

let raises_singular what f =
  if not (try f (); false with L.Singular _ -> true) then fail what

let raises_indefinite what f =
  if not (try f (); false with L.Not_positive_definite _ -> true) then
    fail what

let s_solve () =
  let a = M.of_array 3 3 (mat [| [| 2.0; 1.0; 0.0 |]; [| 1.0; 3.0; 1.0 |]; [| 0.0; 1.0; 2.0 |] |]) in
  let b = M.of_array 3 2 [| 4.0; 10.0; 8.0; 3.0; 5.0; 3.0 |] in
  L.solve a b;
  same "s solve two right hand sides" (M.data b)
    [| 1.0; 2.0; 3.0; 1.0; 1.0; 1.0 |]

let s_lu_solve () =
  let a = M.of_array 3 3 (mat [| [| 2.0; 1.0; 1.0 |]; [| 4.0; -6.0; 0.0 |]; [| -2.0; 7.0; 2.0 |] |]) in
  let ipvt = L.lu a in
  let b = M.of_vec [| 7.0; -8.0; 18.0 |] in
  L.lu_solve a ipvt b;
  same "s lu_solve forward" (M.data b) [| 1.0; 2.0; 3.0 |];
  let c = M.of_vec [| 4.0; 10.0; 7.0 |] in
  L.lu_solve ~trans:true a ipvt c;
  same "s lu_solve transposed" (M.data c) [| 1.0; 2.0; 3.0 |]

let s_lu_rcond () =
  let raw = mat (Array.init 3 (fun i -> Array.init 3 (fun j -> 1.0 /. float_of_int (i + j + 1)))) in
  let a = M.of_array 3 3 (Array.copy raw) in
  let ipvt = Array.make 3 0 in
  let z = Array.make 3 0.0 in
  let r0 = dgeco raw 3 3 ipvt z in
  let got, r1 = L.lu_rcond a in
  near "s lu_rcond estimate" r1 r0;
  same "s lu_rcond factors" (M.data a) raw;
  check "s lu_rcond pivots" (got = ipvt);
  let b = M.of_vec (Array.init 3 (fun i -> (1.0 /. float_of_int (i + 1)) +. (1.0 /. float_of_int (i + 2)) +. (1.0 /. float_of_int (i + 3)))) in
  L.lu_solve a got b;
  same_t 1.0e-8 "s lu_rcond then solve" (M.data b) (ones 3);
  let singular = M.of_array 2 2 (mat [| [| 1.0; 2.0 |]; [| 2.0; 4.0 |] |]) in
  let _, rs = L.lu_rcond singular in
  check "s lu_rcond singular is not an exception" (rs >= 0.0 && rs < 1.0e-14)

let s_band_solve () =
  let abd = M.of_array 4 5 (tridiag_band 5 1 1 4.0 (-1.0)) in
  let b = M.of_array 5 2 [| 3.0; 2.0; 2.0; 2.0; 3.0; 2.0; 4.0; 6.0; 8.0; 16.0 |] in
  L.band_solve ~ml:1 ~mu:1 abd b;
  same "s band_solve two right hand sides" (M.data b)
    [| 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 2.0; 3.0; 4.0; 5.0 |]

let s_band_lu_solve () =
  let abd = M.of_array 4 5 (tridiag_band 5 1 1 4.0 (-1.0)) in
  let raw = tridiag_band 5 1 1 4.0 (-1.0) in
  let ipvt = Array.make 5 0 in
  info_is "s band raw info" (dgbfa raw 4 5 1 1 ipvt) 0;
  let got = L.band_lu ~ml:1 ~mu:1 abd in
  same "s band_lu factors" (M.data abd) raw;
  check "s band_lu pivots" (got = ipvt);
  let b = M.of_vec [| 3.0; 2.0; 2.0; 2.0; 3.0 |] in
  L.band_lu_solve ~ml:1 ~mu:1 abd got b;
  same "s band_lu_solve forward" (M.data b) (ones 5);
  let c = M.of_vec [| 3.0; 2.0; 2.0; 2.0; 3.0 |] in
  L.band_lu_solve ~trans:true ~ml:1 ~mu:1 abd got c;
  same "s band_lu_solve transposed" (M.data c) (ones 5)

let s_cholesky () =
  let rows = [| [| 4.0; 2.0; 0.0 |]; [| 2.0; 5.0; 1.0 |]; [| 0.0; 1.0; 3.0 |] |] in
  let a = M.of_array 3 3 (mat rows) in
  let b = M.of_array 3 2 [| 6.0; 8.0; 4.0; 4.0; 2.0; 0.0 |] in
  L.solve_spd a b;
  same "s solve_spd two right hand sides" (M.data b)
    [| 1.0; 1.0; 1.0; 1.0; 0.0; 0.0 |];
  let f = M.of_array 3 3 (mat rows) in
  L.cholesky f;
  let c = M.of_vec [| 6.0; 8.0; 4.0 |] in
  L.cholesky_solve f c;
  same "s cholesky_solve" (M.data c) (ones 3)

let s_band_cholesky () =
  let abd = M.of_array 2 5 (spd_band 5 1 4.0 (-1.0)) in
  let b = M.of_vec [| 3.0; 2.0; 2.0; 2.0; 3.0 |] in
  L.band_solve_spd ~mu:1 abd b;
  same "s band_solve_spd" (M.data b) (ones 5);
  let f = M.of_array 2 5 (spd_band 5 1 4.0 (-1.0)) in
  L.band_cholesky ~mu:1 f;
  let c = M.of_vec [| 3.0; 2.0; 2.0; 2.0; 3.0 |] in
  L.band_cholesky_solve ~mu:1 f c;
  same "s band_cholesky_solve" (M.data c) (ones 5)

let s_exceptions () =
  raises_singular "s solve singular" (fun () ->
      L.solve
        (M.of_array 2 2 (mat [| [| 1.0; 2.0 |]; [| 2.0; 4.0 |] |]))
        (M.of_vec [| 1.0; 2.0 |]));
  raises_singular "s lu singular" (fun () ->
      ignore (L.lu (M.of_array 2 2 (mat [| [| 1.0; 2.0 |]; [| 2.0; 4.0 |] |]))));
  raises_singular "s band_lu singular" (fun () ->
      ignore (L.band_lu ~ml:1 ~mu:1 (M.create 4 3)));
  raises_indefinite "s cholesky not definite" (fun () ->
      L.cholesky (M.of_array 2 2 (mat [| [| 1.0; 2.0 |]; [| 2.0; 1.0 |] |])));
  raises_indefinite "s solve_spd not definite" (fun () ->
      L.solve_spd
        (M.of_array 2 2 (mat [| [| 1.0; 2.0 |]; [| 2.0; 1.0 |] |]))
        (M.of_vec [| 1.0; 1.0 |]));
  raises_indefinite "s band_cholesky not definite" (fun () ->
      L.band_cholesky ~mu:1 (M.create 2 3))

let s_conformance () =
  raises "s solve not square" (fun () ->
      L.solve (M.create 2 3) (M.of_vec [| 1.0; 2.0 |]));
  raises "s solve rhs rows" (fun () ->
      L.solve (M.create 2 2) (M.of_vec [| 1.0; 2.0; 3.0 |]));
  raises "s cholesky not square" (fun () -> L.cholesky (M.create 2 3));
  raises "s lu_solve short pivots" (fun () ->
      L.lu_solve (M.create 2 2) [| 0 |] (M.of_vec [| 1.0; 2.0 |]));
  raises "s lu_solve pivot out of range" (fun () ->
      L.lu_solve (M.create 2 2) [| 0; 5 |] (M.of_vec [| 1.0; 2.0 |]));
  raises "s band_lu negative bandwidth" (fun () ->
      ignore (L.band_lu ~ml:(-1) ~mu:1 (M.create 4 5)));
  raises "s band_lu too few rows" (fun () ->
      ignore (L.band_lu ~ml:1 ~mu:1 (M.create 3 5)));
  raises "s band_lu bandwidth too wide" (fun () ->
      ignore (L.band_lu ~ml:0 ~mu:3 (M.create 4 3)));
  raises "s band_solve rhs rows" (fun () ->
      L.band_solve ~ml:1 ~mu:1 (M.create 4 5) (M.of_vec [| 1.0; 2.0 |]));
  raises "s band_cholesky too few rows" (fun () ->
      L.band_cholesky ~mu:1 (M.create 1 5))

let () =
  l1_ge_2x2 ();
  l1_ge_identity ();
  l1_ge_hilbert ();
  l1_po_2x2 ();
  l1_po_identity ();
  l1_po_diagonal ();
  l3_ge_doc_example ();
  l3_ge_pascal ();
  l3_ge_vandermonde ();
  l3_po_lehmer ();
  l3_po_tridiagonal ();
  l3_po_identity_plus_ones ();
  l3_ge_forsythe_moler ();
  l3_ge_wilkinson ();
  l3_ge_upper_triangular ();
  l1_gb_tridiagonal ();
  l1_gb_pentadiagonal ();
  l1_gb_diagonal ();
  l1_pb_tridiagonal ();
  l1_pb_diagonal ();
  l1_pb_bandwidth_two ();
  l1_geco_identity ();
  l1_geco_diagonal ();
  l1_geco_hilbert ();
  l2_gb_residual ();
  l2_gb_linearity ();
  l2_gb_inverse_column ();
  l2_gb_zero_rhs ();
  l2_pb_rtr_diagonal ();
  l2_pb_residual ();
  l2_pb_inverse_symmetry ();
  l2_pb_not_definite ();
  l2_geco_range ();
  l2_geco_scaling ();
  l2_geco_near_singular ();
  l2_geco_orthogonal ();
  l2_geco_triangular ();
  l3_gb_laplacian ();
  l3_gb_convection_diffusion ();
  l3_gb_laplacian_slice ();
  l3_pb_fem_bar ();
  l3_pb_fem_mass ();
  l3_pb_string_mode ();
  l3_geco_hilbert ();
  l3_geco_vandermonde ();
  l3_geco_moler ();
  l3_geco_pascal ();
  s_solve ();
  s_lu_solve ();
  s_lu_rcond ();
  s_band_solve ();
  s_band_lu_solve ();
  s_cholesky ();
  s_band_cholesky ();
  s_exceptions ();
  s_conformance ();
  if !bad = 0 then print_string "linpack: all checks passed\n" else exit 1
