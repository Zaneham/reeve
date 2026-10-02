(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Checks for the LAPACK linear system routines. Factorisations are judged by
   rebuilding the original matrix from the factors and solves by the backward
   error of the computed solution, so no reference numbers are needed except
   where the arithmetic is exact. *)

open Reeve.Blasmat
open Reeve.Lapack.Raw

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

let raises_named what name f =
  let m = try f (); "" with Invalid_argument m -> m in
  let k = String.length name in
  check what (String.length m > k && String.sub m 0 k = name)

let seed = ref 1

let rnd () =
  seed := ((!seed * 1103515245) + 12345) mod 2147483648;
  float_of_int ((!seed mod 2001) - 1000) /. 100.0

let gen m n ao lda =
  let a = Array.make (ao + (lda * n)) 0.0 in
  for j = 0 to n - 1 do
    for i = 0 to m - 1 do
      a.(ao + i + (j * lda)) <- rnd ()
    done
  done;
  a

let spd n ao lda =
  let a = Array.make (ao + (lda * n)) 0.0 in
  for j = 0 to n - 1 do
    for i = j + 1 to n - 1 do
      let v = rnd () in
      a.(ao + i + (j * lda)) <- v;
      a.(ao + j + (i * lda)) <- v
    done
  done;
  for j = 0 to n - 1 do
    a.(ao + j + (j * lda)) <- (10.0 *. float_of_int n) +. float_of_int (j + 1)
  done;
  a

let amax a ao lda m n =
  let v = ref 0.0 in
  for j = 0 to n - 1 do
    for i = 0 to m - 1 do
      v := Float.max !v (Float.abs a.(ao + i + (j * lda)))
    done
  done;
  !v

let residual what trans a ao lda n x xo ldx b bo ldb nrhs =
  let an = amax a ao lda n n in
  let xn = amax x xo ldx n nrhs in
  let worst = ref 0.0 in
  for j = 0 to nrhs - 1 do
    for i = 0 to n - 1 do
      let s = ref 0.0 in
      for k = 0 to n - 1 do
        let aik =
          match trans with
          | No_trans -> a.(ao + i + (k * lda))
          | Trans | Conj_trans -> a.(ao + k + (i * lda))
        in
        s := !s +. (aik *. x.(xo + k + (j * ldx)))
      done;
      worst := Float.max !worst (Float.abs (!s -. b.(bo + i + (j * ldb))))
    done
  done;
  let bound = 1.0e-11 *. float_of_int (n + 1) *. (1.0 +. an) *. (1.0 +. xn) in
  check (Printf.sprintf "%s residual" what) (!worst < bound)

let lu_rebuild a ao lda m n ipiv =
  let mn = min m n in
  let r = Array.make (lda * n) 0.0 in
  for j = 0 to n - 1 do
    for i = 0 to m - 1 do
      let s = ref 0.0 in
      for k = 0 to mn - 1 do
        let lik =
          if k < i then a.(ao + i + (k * lda))
          else if k = i then 1.0
          else 0.0
        in
        let ukj = if k <= j then a.(ao + k + (j * lda)) else 0.0 in
        s := !s +. (lik *. ukj)
      done;
      r.(i + (j * lda)) <- !s
    done
  done;
  for k = mn - 1 downto 0 do
    let p = ipiv.(k) in
    if p <> k then
      for j = 0 to n - 1 do
        let t = r.(k + (j * lda)) in
        r.(k + (j * lda)) <- r.(p + (j * lda));
        r.(p + (j * lda)) <- t
      done
  done;
  r

let lu_rebuild_ok what orig ao lda m n a ipiv =
  let r = lu_rebuild a ao lda m n ipiv in
  let an = amax orig ao lda m n in
  let worst = ref 0.0 in
  for j = 0 to n - 1 do
    for i = 0 to m - 1 do
      worst :=
        Float.max !worst
          (Float.abs (r.(i + (j * lda)) -. orig.(ao + i + (j * lda))))
    done
  done;
  let bound = 1.0e-11 *. float_of_int (min m n + 1) *. (1.0 +. an) in
  check (Printf.sprintf "%s rebuilds a" what) (!worst < bound)

let chol_rebuild_ok what uplo orig a ao lda n =
  let an = amax orig ao lda n n in
  let worst = ref 0.0 in
  for j = 0 to n - 1 do
    let lo = match uplo with Upper -> 0 | Lower -> j in
    let hi = match uplo with Upper -> j | Lower -> n - 1 in
    for i = lo to hi do
      let s = ref 0.0 in
      for k = 0 to min i j do
        let p =
          match uplo with
          | Upper -> a.(ao + k + (i * lda)) *. a.(ao + k + (j * lda))
          | Lower -> a.(ao + i + (k * lda)) *. a.(ao + j + (k * lda))
        in
        s := !s +. p
      done;
      worst := Float.max !worst (Float.abs (!s -. orig.(ao + i + (j * lda))))
    done
  done;
  let bound = 1.0e-11 *. float_of_int (n + 1) *. (1.0 +. an) in
  check (Printf.sprintf "%s rebuilds a" what) (!worst < bound)

let sizes = [ 1; 2; 3; 4; 5; 8; 13; 31; 63; 64; 65; 80; 100 ]
let shapes = [ (1, 1); (1, 4); (4, 1); (3, 5); (5, 3); (8, 8); (13, 7); (7, 13); (65, 40); (40, 65); (70, 70) ]

let gesv_diagonal_exact () =
  let a = mat [| [| 2.0; 0.0; 0.0 |]; [| 0.0; 4.0; 0.0 |]; [| 0.0; 0.0; 8.0 |] |] in
  let b = [| 2.0; 8.0; 24.0 |] in
  let ipiv = Array.make 3 0 in
  info_is "gesv diagonal info" (dgesv 3 1 a 0 3 ipiv 0 b 0 3) 0;
  same "gesv diagonal" b [| 1.0; 2.0; 3.0 |]

let gesv_small_exact () =
  let a = mat [| [| 1.0; 2.0 |]; [| 3.0; 4.0 |] |] in
  let b = [| 5.0; 11.0 |] in
  let ipiv = Array.make 2 0 in
  info_is "gesv small info" (dgesv 2 1 a 0 2 ipiv 0 b 0 2) 0;
  same "gesv small" b [| 1.0; 2.0 |]

let gesv_two_right_hand_sides () =
  let a = mat [| [| 2.0; 0.0 |]; [| 0.0; 4.0 |] |] in
  let b = mat [| [| 2.0; 6.0 |]; [| 8.0; 12.0 |] |] in
  let ipiv = Array.make 2 0 in
  info_is "gesv nrhs info" (dgesv 2 2 a 0 2 ipiv 0 b 0 2) 0;
  same "gesv nrhs" b (mat [| [| 1.0; 3.0 |]; [| 2.0; 3.0 |] |])

let gesv_residual () =
  seed := 3;
  List.iter
    (fun n ->
      List.iter
        (fun nrhs ->
          let lda = n + 2 in
          let a = gen n n 0 lda in
          let orig = Array.copy a in
          let b = gen n nrhs 0 lda in
          let rhs = Array.copy b in
          let ipiv = Array.make n 0 in
          let info = dgesv n nrhs a 0 lda ipiv 0 b 0 lda in
          info_is (Printf.sprintf "gesv %d x %d info" n nrhs) info 0;
          residual
            (Printf.sprintf "gesv %d x %d" n nrhs)
            No_trans orig 0 lda n b 0 lda rhs 0 lda nrhs)
        [ 1; 3 ])
    sizes

let gesv_offset_leaves_padding () =
  seed := 5;
  let n = 6 and lda = 9 and ao = 4 in
  let a = gen n n ao lda in
  Array.iteri (fun i v -> if v = 0.0 then a.(i) <- 77.0 +. float_of_int i) a;
  let orig = Array.copy a in
  let b = gen n 2 ao lda in
  Array.iteri (fun i v -> if v = 0.0 then b.(i) <- 88.0 +. float_of_int i) b;
  let rhs = Array.copy b in
  let ipiv = Array.make n 0 in
  info_is "gesv offset info" (dgesv n 2 a ao lda ipiv 0 b ao lda) 0;
  residual "gesv offset" No_trans orig ao lda n b ao lda rhs ao lda 2;
  let touched = ref 0 in
  for i = 0 to ao - 1 do
    if b.(i) <> rhs.(i) then touched := !touched + 1
  done;
  for j = 0 to 1 do
    for i = n to lda - 1 do
      if b.(ao + i + (j * lda)) <> rhs.(ao + i + (j * lda)) then
        touched := !touched + 1
    done
  done;
  check "gesv offset leaves padding" (!touched = 0)

let gesv_singular () =
  let a = mat [| [| 1.0; 2.0 |]; [| 2.0; 4.0 |] |] in
  let b = [| 1.0; 1.0 |] in
  let keep = Array.copy b in
  let ipiv = Array.make 2 0 in
  info_is "gesv singular info" (dgesv 2 1 a 0 2 ipiv 0 b 0 2) 2;
  same "gesv singular leaves b" b keep;
  let a = mat [| [| 0.0; 0.0 |]; [| 0.0; 1.0 |] |] in
  let ipiv = Array.make 2 0 in
  info_is "gesv zero column info" (dgesv 2 1 a 0 2 ipiv 0 b 0 2) 1;
  let a =
    mat
      [| [| 1.0; 0.0; 2.0 |]; [| 3.0; 0.0; 4.0 |]; [| 5.0; 0.0; 6.0 |] |]
  in
  let b3 = [| 1.0; 1.0; 1.0 |] in
  let ipiv = Array.make 3 0 in
  info_is "gesv middle zero column info" (dgesv 3 1 a 0 3 ipiv 0 b3 0 3) 2

let getf2_rebuild () =
  seed := 13;
  List.iter
    (fun (m, n) ->
      let lda = m + 1 in
      let a = gen m n 0 lda in
      let orig = Array.copy a in
      let ipiv = Array.make (max (min m n) 1) 0 in
      let info = dgetf2 m n a 0 lda ipiv 0 in
      info_is (Printf.sprintf "getf2 %d x %d info" m n) info 0;
      lu_rebuild_ok (Printf.sprintf "getf2 %d x %d" m n) orig 0 lda m n a ipiv)
    shapes

let getrf2_rebuild () =
  seed := 17;
  List.iter
    (fun (m, n) ->
      let lda = m + 1 in
      let a = gen m n 0 lda in
      let orig = Array.copy a in
      let ipiv = Array.make (max (min m n) 1) 0 in
      let info = dgetrf2 m n a 0 lda ipiv 0 in
      info_is (Printf.sprintf "getrf2 %d x %d info" m n) info 0;
      lu_rebuild_ok (Printf.sprintf "getrf2 %d x %d" m n) orig 0 lda m n a ipiv)
    shapes

let getrf_rebuild () =
  seed := 19;
  List.iter
    (fun (m, n) ->
      let lda = m + 3 in
      let a = gen m n 0 lda in
      let orig = Array.copy a in
      let ipiv = Array.make (max (min m n) 1) 0 in
      let info = dgetrf m n a 0 lda ipiv 0 in
      info_is (Printf.sprintf "getrf %d x %d info" m n) info 0;
      lu_rebuild_ok (Printf.sprintf "getrf %d x %d" m n) orig 0 lda m n a ipiv)
    shapes

let getf2_singular_info () =
  let a = mat [| [| 1.0; 2.0 |]; [| 2.0; 4.0 |] |] in
  let ipiv = Array.make 2 0 in
  info_is "getf2 singular info" (dgetf2 2 2 a 0 2 ipiv 0) 2;
  let a = mat [| [| 0.0; 0.0 |]; [| 0.0; 1.0 |] |] in
  info_is "getf2 zero column info" (dgetf2 2 2 a 0 2 ipiv 0) 1;
  let a = mat [| [| 1.0; 2.0 |]; [| 2.0; 4.0 |] |] in
  info_is "getrf2 singular info" (dgetrf2 2 2 a 0 2 ipiv 0) 2;
  let a = mat [| [| 0.0; 0.0 |]; [| 0.0; 1.0 |] |] in
  info_is "getrf2 zero column info" (dgetrf2 2 2 a 0 2 ipiv 0) 1;
  let a = mat [| [| 1.0; 2.0 |]; [| 2.0; 4.0 |] |] in
  info_is "getrf singular info" (dgetrf 2 2 a 0 2 ipiv 0) 2

let ipiv_is_zero_based () =
  let a = mat [| [| 0.0; 1.0 |]; [| 1.0; 0.0 |] |] in
  let ipiv = Array.make 2 (-1) in
  info_is "ipiv zero based info" (dgetrf 2 2 a 0 2 ipiv 0) 0;
  check "ipiv zero based first" (ipiv.(0) = 1);
  check "ipiv zero based second" (ipiv.(1) = 1);
  let a = mat [| [| 2.0; 1.0 |]; [| 1.0; 2.0 |] |] in
  let ipiv = Array.make 2 (-1) in
  info_is "ipiv no swap info" (dgetrf 2 2 a 0 2 ipiv 0) 0;
  check "ipiv no swap first" (ipiv.(0) = 0);
  check "ipiv no swap second" (ipiv.(1) = 1)

let ipiv_offset () =
  let a = mat [| [| 0.0; 1.0 |]; [| 1.0; 0.0 |] |] in
  let ipiv = Array.make 5 (-1) in
  info_is "ipiv offset info" (dgetrf 2 2 a 0 2 ipiv 3) 0;
  check "ipiv offset untouched" (ipiv.(0) = -1 && ipiv.(1) = -1 && ipiv.(2) = -1);
  check "ipiv offset written" (ipiv.(3) = 1 && ipiv.(4) = 1)

let getrs_both_ways () =
  seed := 23;
  List.iter
    (fun n ->
      List.iter
        (fun t ->
          let lda = n + 1 in
          let a = gen n n 0 lda in
          let orig = Array.copy a in
          let b = gen n 2 0 lda in
          let rhs = Array.copy b in
          let ipiv = Array.make n 0 in
          info_is (Printf.sprintf "getrs %d factor info" n)
            (dgetrf n n a 0 lda ipiv 0) 0;
          info_is (Printf.sprintf "getrs %d solve info" n)
            (dgetrs t n 2 a 0 lda ipiv 0 b 0 lda) 0;
          residual
            (Printf.sprintf "getrs %d %s" n
               (match t with
               | No_trans -> "notrans"
               | Trans -> "trans"
               | Conj_trans -> "conjtrans"))
            t orig 0 lda n b 0 lda rhs 0 lda 2)
        [ No_trans; Trans; Conj_trans ])
    [ 1; 2; 3; 5; 8; 13; 65; 80 ]

let getrs_after_getf2 () =
  seed := 29;
  let n = 7 in
  let lda = n + 2 in
  let a = gen n n 0 lda in
  let orig = Array.copy a in
  let b = gen n 3 0 lda in
  let rhs = Array.copy b in
  let ipiv = Array.make n 0 in
  info_is "getf2 then getrs factor info" (dgetf2 n n a 0 lda ipiv 0) 0;
  info_is "getf2 then getrs solve info"
    (dgetrs No_trans n 3 a 0 lda ipiv 0 b 0 lda) 0;
  residual "getf2 then getrs" No_trans orig 0 lda n b 0 lda rhs 0 lda 3

let potf2_diagonal_exact () =
  let a = mat [| [| 4.0; 0.0; 0.0 |]; [| 0.0; 9.0; 0.0 |]; [| 0.0; 0.0; 16.0 |] |] in
  info_is "potf2 diagonal info" (dpotf2 Upper 3 a 0 3) 0;
  near "potf2 diagonal 11" a.(0) 2.0;
  near "potf2 diagonal 22" a.(4) 3.0;
  near "potf2 diagonal 33" a.(8) 4.0;
  let b = [| 8.0; 27.0; 64.0 |] in
  info_is "potrs diagonal info" (dpotrs Upper 3 1 a 0 3 b 0 3) 0;
  same "potrs diagonal" b [| 2.0; 3.0; 4.0 |]

let posv_diagonal_exact () =
  let a = mat [| [| 4.0; 0.0 |]; [| 0.0; 16.0 |] |] in
  let b = [| 8.0; 64.0 |] in
  info_is "posv diagonal info" (dposv Lower 2 1 a 0 2 b 0 2) 0;
  same "posv diagonal" b [| 2.0; 4.0 |]

let chol_rebuild which () =
  seed := 31 + which;
  List.iter
    (fun n ->
      List.iter
        (fun u ->
          let lda = n + 2 in
          let a = spd n 0 lda in
          let orig = Array.copy a in
          let label =
            Printf.sprintf "%s %d %s"
              (if which = 0 then "potf2"
               else if which = 1 then "potrf2"
               else "potrf")
              n
              (match u with Upper -> "upper" | Lower -> "lower")
          in
          let info =
            if which = 0 then dpotf2 u n a 0 lda
            else if which = 1 then dpotrf2 u n a 0 lda
            else dpotrf u n a 0 lda
          in
          info_is (label ^ " info") info 0;
          chol_rebuild_ok label u orig a 0 lda n)
        [ Upper; Lower ])
    sizes

let posv_residual () =
  seed := 41;
  List.iter
    (fun n ->
      List.iter
        (fun nrhs ->
          List.iter
            (fun u ->
              let lda = n + 1 in
              let a = spd n 0 lda in
              let orig = Array.copy a in
              let b = gen n nrhs 0 lda in
              let rhs = Array.copy b in
              let label =
                Printf.sprintf "posv %d x %d %s" n nrhs
                  (match u with Upper -> "upper" | Lower -> "lower")
              in
              info_is (label ^ " info") (dposv u n nrhs a 0 lda b 0 lda) 0;
              residual label No_trans orig 0 lda n b 0 lda rhs 0 lda nrhs)
            [ Upper; Lower ])
        [ 1; 3 ])
    sizes

let potrs_after_potf2 () =
  seed := 47;
  let n = 9 in
  let lda = n + 2 in
  List.iter
    (fun u ->
      let a = spd n 0 lda in
      let orig = Array.copy a in
      let b = gen n 2 0 lda in
      let rhs = Array.copy b in
      info_is "potf2 then potrs factor info" (dpotf2 u n a 0 lda) 0;
      info_is "potf2 then potrs solve info" (dpotrs u n 2 a 0 lda b 0 lda) 0;
      residual "potf2 then potrs" No_trans orig 0 lda n b 0 lda rhs 0 lda 2)
    [ Upper; Lower ]

let chol_not_definite () =
  List.iter
    (fun u ->
      let a = mat [| [| 1.0; 2.0 |]; [| 2.0; 1.0 |] |] in
      info_is "potf2 indefinite info" (dpotf2 u 2 a 0 2) 2;
      let a = mat [| [| 1.0; 2.0 |]; [| 2.0; 1.0 |] |] in
      info_is "potrf2 indefinite info" (dpotrf2 u 2 a 0 2) 2;
      let a = mat [| [| 1.0; 2.0 |]; [| 2.0; 1.0 |] |] in
      info_is "potrf indefinite info" (dpotrf u 2 a 0 2) 2;
      let a = mat [| [| -1.0 |] |] in
      info_is "potf2 negative info" (dpotf2 u 1 a 0 1) 1;
      let a = mat [| [| 0.0 |] |] in
      info_is "potrf2 zero info" (dpotrf2 u 1 a 0 1) 1;
      let a =
        mat
          [| [| 4.0; 0.0; 0.0 |]; [| 0.0; 4.0; 0.0 |]; [| 0.0; 0.0; -1.0 |] |]
      in
      info_is "potrf third minor info" (dpotrf u 3 a 0 3) 3;
      let a = mat [| [| 4.0; 0.0 |]; [| 0.0; -1.0 |] |] in
      let b = [| 1.0; 1.0 |] in
      let keep = Array.copy b in
      info_is "posv indefinite info" (dposv u 2 1 a 0 2 b 0 2) 2;
      same "posv indefinite leaves b" b keep)
    [ Upper; Lower ]

let chol_blocked_boundary () =
  seed := 53;
  List.iter
    (fun n ->
      List.iter
        (fun u ->
          let lda = n in
          let a = spd n 0 lda in
          let orig = Array.copy a in
          info_is
            (Printf.sprintf "potrf boundary %d info" n)
            (dpotrf u n a 0 lda) 0;
          chol_rebuild_ok
            (Printf.sprintf "potrf boundary %d" n)
            u orig a 0 lda n)
        [ Upper; Lower ])
    [ 62; 63; 64; 65; 66; 127; 128; 129 ]

let lu_blocked_boundary () =
  seed := 59;
  List.iter
    (fun n ->
      let lda = n in
      let a = gen n n 0 lda in
      let orig = Array.copy a in
      let ipiv = Array.make n 0 in
      info_is
        (Printf.sprintf "getrf boundary %d info" n)
        (dgetrf n n a 0 lda ipiv 0) 0;
      lu_rebuild_ok
        (Printf.sprintf "getrf boundary %d" n)
        orig 0 lda n n a ipiv)
    [ 62; 63; 64; 65; 66; 127; 128; 129 ]

let zero_sizes () =
  let a = [| 1.0 |] and b = [| 2.0 |] in
  let ipiv = Array.make 1 (-1) in
  info_is "gesv n zero" (dgesv 0 1 a 0 1 ipiv 0 b 0 1) 0;
  info_is "gesv nrhs zero" (dgesv 1 0 a 0 1 ipiv 0 b 0 1) 0;
  info_is "getrf m zero" (dgetrf 0 1 a 0 1 ipiv 0) 0;
  info_is "getrf n zero" (dgetrf 1 0 a 0 1 ipiv 0) 0;
  info_is "getrf2 m zero" (dgetrf2 0 1 a 0 1 ipiv 0) 0;
  info_is "getf2 n zero" (dgetf2 1 0 a 0 1 ipiv 0) 0;
  info_is "getrs n zero" (dgetrs No_trans 0 1 a 0 1 ipiv 0 b 0 1) 0;
  info_is "getrs nrhs zero" (dgetrs No_trans 1 0 a 0 1 ipiv 0 b 0 1) 0;
  info_is "potrf n zero" (dpotrf Upper 0 a 0 1) 0;
  info_is "potrf2 n zero" (dpotrf2 Upper 0 a 0 1) 0;
  info_is "potf2 n zero" (dpotf2 Upper 0 a 0 1) 0;
  info_is "potrs n zero" (dpotrs Upper 0 1 a 0 1 b 0 1) 0;
  info_is "potrs nrhs zero" (dpotrs Upper 1 0 a 0 1 b 0 1) 0;
  info_is "posv n zero" (dposv Upper 0 1 a 0 1 b 0 1) 0;
  same "zero sizes leave a" a [| 1.0 |];
  same "zero sizes leave b" b [| 2.0 |]

let bad_arguments () =
  let one = [| 1.0 |] in
  let ip = Array.make 1 0 in
  raises_named "gesv bad n" "dgesv" (fun () ->
      ignore (dgesv (-1) 1 one 0 1 ip 0 one 0 1));
  raises_named "gesv bad nrhs" "dgesv" (fun () ->
      ignore (dgesv 1 (-1) one 0 1 ip 0 one 0 1));
  raises_named "gesv bad lda" "dgesv" (fun () ->
      ignore (dgesv 2 1 one 0 1 ip 0 one 0 2));
  raises_named "gesv bad ldb" "dgesv" (fun () ->
      ignore (dgesv 2 1 (Array.make 4 0.0) 0 2 ip 0 one 0 1));
  raises_named "getrf bad m" "dgetrf" (fun () ->
      ignore (dgetrf (-1) 1 one 0 1 ip 0));
  raises_named "getrf bad lda" "dgetrf" (fun () ->
      ignore (dgetrf 2 1 one 0 1 ip 0));
  raises_named "getrf2 bad n" "dgetrf2" (fun () ->
      ignore (dgetrf2 1 (-1) one 0 1 ip 0));
  raises_named "getf2 bad lda" "dgetf2" (fun () ->
      ignore (dgetf2 2 1 one 0 1 ip 0));
  raises_named "getrs bad n" "dgetrs" (fun () ->
      ignore (dgetrs No_trans (-1) 1 one 0 1 ip 0 one 0 1));
  raises_named "getrs bad ldb" "dgetrs" (fun () ->
      ignore (dgetrs No_trans 2 1 (Array.make 4 0.0) 0 2 ip 0 one 0 1));
  raises_named "potrf bad n" "dpotrf" (fun () ->
      ignore (dpotrf Upper (-1) one 0 1));
  raises_named "potrf bad lda" "dpotrf" (fun () ->
      ignore (dpotrf Upper 2 one 0 1));
  raises_named "potrf2 bad n" "dpotrf2" (fun () ->
      ignore (dpotrf2 Lower (-1) one 0 1));
  raises_named "potf2 bad lda" "dpotf2" (fun () ->
      ignore (dpotf2 Lower 2 one 0 1));
  raises_named "potrs bad nrhs" "dpotrs" (fun () ->
      ignore (dpotrs Upper 1 (-1) one 0 1 one 0 1));
  raises_named "potrs bad ldb" "dpotrs" (fun () ->
      ignore (dpotrs Upper 2 1 (Array.make 4 0.0) 0 2 one 0 1));
  raises_named "posv bad n" "dposv" (fun () ->
      ignore (dposv Upper (-1) 1 one 0 1 one 0 1));
  raises_named "posv bad ldb" "dposv" (fun () ->
      ignore (dposv Upper 2 1 (Array.make 4 0.0) 0 2 one 0 1))

let () =
  gesv_diagonal_exact ();
  gesv_small_exact ();
  gesv_two_right_hand_sides ();
  gesv_residual ();
  gesv_offset_leaves_padding ();
  gesv_singular ();
  getf2_rebuild ();
  getrf2_rebuild ();
  getrf_rebuild ();
  getf2_singular_info ();
  ipiv_is_zero_based ();
  ipiv_offset ();
  getrs_both_ways ();
  getrs_after_getf2 ();
  potf2_diagonal_exact ();
  posv_diagonal_exact ();
  chol_rebuild 0 ();
  chol_rebuild 1 ();
  chol_rebuild 2 ();
  posv_residual ();
  potrs_after_potf2 ();
  chol_not_definite ();
  chol_blocked_boundary ();
  lu_blocked_boundary ();
  zero_sizes ();
  bad_arguments ();
  if !bad = 0 then print_string "lapack: all checks passed\n" else exit 1
