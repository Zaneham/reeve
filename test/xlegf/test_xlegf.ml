(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Level 1, level 2 and level 3 checks for the extended-range arithmetic and
   the Legendre functions of the first and second kind built on it. *)

open Reeve.Xlegf

let pi = 3.14159265358979323846
let env = dxset 0 0 0.0 0
let tol = 1.0e-13
let bad = ref 0

let fail what =
  Printf.printf "FAIL %s\n" what;
  bad := !bad + 1

let holds what c = if not c then fail what
let near_t what t a b = if not (Float.abs (a -. b) < t) then fail what

let rel_t what t a b =
  let e =
    if Float.abs b > 1.0e-14 then Float.abs (a -. b) /. Float.abs b
    else Float.abs (a -. b)
  in
  if not (e < t) then fail what

let raises what f =
  match f () with
  | () -> fail what
  | exception Invalid_argument _ -> ()

let same what a b = holds what (a.x = b.x && a.ix = b.ix)
let log2v a = Float.log2 (Float.abs a.x) +. float_of_int a.ix
let sgn m = if m land 1 = 0 then 1.0 else -1.0

let fact n =
  let r = ref 1.0 in
  for i = 2 to n do
    r := !r *. float_of_int i
  done;
  !r

let dfact n =
  let r = ref 1.0 and i = ref n in
  while !i > 1 do
    r := !r *. float_of_int !i;
    i := !i - 2
  done;
  !r

let pnm0 n m =
  if (n - m) mod 2 <> 0 then 0.0
  else sgn ((n - m) / 2) *. dfact (n + m - 1) /. dfact (n - m)

let legp n x =
  match n with
  | 0 -> 1.0
  | 1 -> x
  | 2 -> ((3.0 *. x *. x) -. 1.0) /. 2.0
  | 3 -> ((5.0 *. x *. x *. x) -. (3.0 *. x)) /. 2.0
  | _ ->
    ((35.0 *. x *. x *. x *. x) -. (30.0 *. x *. x) +. 3.0) /. 8.0

let atanh2 x = 0.5 *. log ((1.0 +. x) /. (1.0 -. x))

(* ---- Level 1, the arithmetic ---- *)

let normals =
  [| 1.0; -1.0; 0.5; -3.5; 123.456; 1.0e-300; 1.0e300; min_float; 0.1 |]

let l1_roundtrip () =
  Array.iter
    (fun v ->
      let a = dxadj env { x = v; ix = 0 } in
      holds "l1 adjusted principal part in range"
        (Float.abs a.x >= 1.0 /. 6.7e153 && Float.abs a.x < 7.0e153);
      near_t "l1 dxadj keeps the value" 1.0e-9 (log2v a) (Float.log2 (Float.abs v));
      let b = dxred env a in
      holds "l1 dxred recovers the float" (b.x = v && b.ix = 0);
      let c = dxred env { x = v; ix = 0 } in
      holds "l1 dxred of a reduced number" (c.x = v && c.ix = 0))
    normals

let l1_zero () =
  same "l1 dxadj zero" (dxadj env { x = 0.0; ix = 77 }) { x = 0.0; ix = 0 };
  same "l1 dxred zero" (dxred env { x = 0.0; ix = 77 }) { x = 0.0; ix = 0 };
  same "l1 dxadd zero left" (dxadd env { x = 0.0; ix = 5 } { x = 3.0; ix = 0 })
    { x = 3.0; ix = 0 };
  same "l1 dxadd zero right" (dxadd env { x = 3.0; ix = 0 } { x = 0.0; ix = 5 })
    { x = 3.0; ix = 0 }

let add_pairs =
  [|
    (1.0, 2.0);
    (1.5, -0.25);
    (1.0e200, 1.0e-200);
    (-3.0, 7.0);
    (0.1, 0.2);
    (1.0e-300, 1.0e-300);
    (1.0e308, -1.0e308);
    (123.456, -0.000789);
    (1.0e300, 1.0e300);
  |]

let l1_add_matches_float () =
  Array.iter
    (fun (a, b) ->
      let s = dxred env (dxadd env { x = a; ix = 0 } { x = b; ix = 0 }) in
      holds "l1 dxadd agrees with float addition" (s.x = a +. b && s.ix = 0))
    add_pairs

let l1_extended () =
  same "l1 dxadd shares an index"
    (dxadd env { x = 1.0; ix = 3000 } { x = 1.0; ix = 2998 })
    { x = 1.25; ix = 3000 };
  same "l1 dxadd drops a negligible operand"
    (dxadd env { x = 1.0; ix = 4000 } { x = 1.0; ix = -4000 })
    { x = 1.0; ix = 4000 };
  same "l1 dxred leaves an unrepresentable value alone"
    (dxred env { x = 1.0; ix = -2000 })
    { x = 1.0; ix = -2000 };
  let tiny = dxadj env { x = 1.0e-320; ix = 0 } in
  let back = dxred env tiny in
  holds "l1 dxred leaves a subnormal extended" (back.ix <> 0);
  near_t "l1 dxadj keeps a subnormal value" 1.0e-9 (log2v tiny)
    (Float.log2 1.0e-320)

let con_indices = [| 1; -1; 100; -100; 5000; -5000; 123457; -123457 |]

let l1_dxcon () =
  Array.iter
    (fun k ->
      let a = dxcon env { x = 1.0; ix = k } in
      rel_t "l1 dxcon preserves the value" 1.0e-14
        (log10 (Float.abs a.x) +. float_of_int a.ix)
        (float_of_int k *. log10 2.0);
      if a.ix <> 0 then
        holds "l1 dxcon normalises to a decimal fraction"
          (Float.abs a.x >= 0.1 && Float.abs a.x < 1.0))
    con_indices;
  same "l1 dxcon leaves a reduced number alone"
    (dxcon env { x = 0.625; ix = 0 })
    { x = 0.625; ix = 0 }

(* ---- Level 1, the functions ---- *)

let l1_legendre_closed () =
  List.iter
    (fun th ->
      let x = cos th in
      let p = dxlegf env 0.0 4 0 0 th Pneg in
      let r = dxlegf env 0.0 4 0 0 th Ppos in
      for n = 0 to 4 do
        rel_t "l1 P of degree n" tol p.(n).x (legp n x);
        holds "l1 P needs no index" (p.(n).ix = 0);
        rel_t "l1 P of positive order with mu zero" tol r.(n).x (legp n x)
      done;
      let q0 = dxlegf env 0.0 0 0 0 th Q in
      let q1 = dxlegf env 1.0 0 0 0 th Q in
      rel_t "l1 Q of degree zero" tol q0.(0).x (atanh2 x);
      rel_t "l1 Q of degree one" tol q1.(0).x ((x *. atanh2 x) -. 1.0))
    [ 1.0; 0.4; 1.2; pi /. 2.0 ]

let l1_legendre_at_zero () =
  let th = pi /. 2.0 in
  for n = 2 to 4 do
    let dn = float_of_int n in
    let pneg = dxlegf env dn 0 0 n th Pneg in
    let ppos = dxlegf env dn 0 0 n th Ppos in
    let pnrm = dxlegf env dn 0 0 n th Pnorm in
    let nrmp, _ = dxnrmp env n 0 n 0.0 By_x in
    for m = 0 to n do
      let p = pnm0 n m in
      let f = fact (n - m) /. fact (n + m) in
      let name s = Printf.sprintf "l1 %s(%d,%d) at x = 0" s n m in
      rel_t (name "Pneg") tol pneg.(m).x (f *. p);
      rel_t (name "Ppos") tol ppos.(m).x (sgn m *. p);
      rel_t (name "Pnorm") tol pnrm.(m).x (sqrt ((dn +. 0.5) *. f) *. p);
      rel_t (name "dxnrmp") tol nrmp.(m).x (sqrt ((dn +. 0.5) *. f) *. p)
    done
  done

let l1_nrmp_ends () =
  let p, isig = dxnrmp env 5 0 0 1.0 By_x in
  rel_t "l1 dxnrmp at x = 1" tol p.(0).x (sqrt 5.5);
  holds "l1 dxnrmp isig at x = 1" (isig = 1);
  let p, _ = dxnrmp env 5 0 0 (-1.0) By_x in
  rel_t "l1 dxnrmp at x = -1, odd degree" tol p.(0).x (-.sqrt 5.5);
  let p, _ = dxnrmp env 4 0 0 (-1.0) By_x in
  rel_t "l1 dxnrmp at x = -1, even degree" tol p.(0).x (sqrt 4.5);
  let p, _ = dxnrmp env 5 0 0 0.0 By_theta in
  rel_t "l1 dxnrmp at theta = 0" tol p.(0).x (sqrt 5.5);
  let p, _ = dxnrmp env 0 0 0 0.3 By_x in
  rel_t "l1 dxnrmp of degree zero" tol p.(0).x (sqrt 0.5);
  let p, isig = dxnrmp env 5 1 2 1.0 By_x in
  holds "l1 dxnrmp at x = 1 above order zero"
    (p.(0).x = 0.0 && p.(1).x = 0.0 && isig = 0);
  let p, isig = dxnrmp env 2 3 4 0.3 By_x in
  holds "l1 dxnrmp vanishes above the degree"
    (p.(0).x = 0.0 && p.(1).x = 0.0 && isig = 0)

(* ---- Level 2, identities ---- *)

let l2_arith () =
  Array.iter
    (fun (a, b) ->
      let x = { x = a; ix = 11 } and y = { x = b; ix = -7 } in
      same "l2 dxadd commutes" (dxadd env x y) (dxadd env y x))
    add_pairs;
  let a = { x = 1.0; ix = 700 } and b = { x = 1.0; ix = 700 } in
  let s = dxadd env a b in
  near_t "l2 dxadd doubles" 1.0e-9 (log2v s) (701.0);
  let c = { x = 3.0; ix = -4000 } in
  let d = dxadd env c { x = -1.0; ix = -4000 } in
  near_t "l2 dxadd subtracts in the extended range" 1.0e-9 (log2v d)
    (Float.log2 2.0 -. 4000.0)

let l2_degree_recurrence () =
  let th = 0.7 in
  let x = cos th in
  List.iter
    (fun kind ->
      let p = dxlegf env 0.0 6 0 0 th kind in
      for n = 1 to 5 do
        let dn = float_of_int n in
        let lhs = (dn +. 1.0) *. p.(n + 1).x in
        let rhs =
          (((2.0 *. dn) +. 1.0) *. x *. p.(n).x) -. (dn *. p.(n - 1).x)
        in
        rel_t "l2 degree recurrence" 1.0e-12 lhs rhs
      done)
    [ Pneg; Q ]

let l2_nrmp_modes () =
  let th = 0.9 in
  let a, sa = dxnrmp env 12 0 12 (cos th) By_x in
  let b, sb = dxnrmp env 12 0 12 th By_theta in
  holds "l2 dxnrmp isig agrees between modes" (sa = sb);
  Array.iteri
    (fun i v ->
      rel_t "l2 dxnrmp agrees between modes" 1.0e-14 v.x b.(i).x;
      holds "l2 dxnrmp index agrees between modes" (v.ix = b.(i).ix))
    a

let l2_nrmp_against_legf () =
  let th = 0.9 in
  for n = 1 to 20 do
    let p, _ = dxnrmp env n 0 0 (cos th) By_x in
    let dn = float_of_int n in
    if n <= 4 then
      rel_t "l2 dxnrmp of order zero" tol p.(0).x
        (sqrt (dn +. 0.5) *. legp n (cos th));
    let q = dxlegf env dn 0 0 0 th Pnorm in
    rel_t "l2 dxnrmp agrees with dxlegf" 1.0e-12 p.(0).x q.(0).x
  done

let mul a b = dxadj env { x = a.x *. b.x; ix = a.ix + b.ix }
let sub a b = dxadd env a { x = -.b.x; ix = b.ix }

let casorati_nu dnu1 nudiff mu1 theta =
  let p = dxlegf env dnu1 nudiff mu1 mu1 theta Pneg in
  let q = dxlegf env dnu1 nudiff mu1 mu1 theta Q in
  let r = dxlegf env dnu1 nudiff mu1 mu1 theta Ppos in
  let dmu1 = float_of_int mu1 in
  for i = 0 to nudiff - 1 do
    let nu = dnu1 +. float_of_int i in
    let x1 = mul p.(i + 1) q.(i) in
    let x2 = mul p.(i) q.(i + 1) in
    let x1 = { x1 with x = (dmu1 +. nu +. 1.0) *. x1.x } in
    let x2 = { x2 with x = (dmu1 -. nu -. 1.0) *. x2.x } in
    let c = dxadj env (dxadd env x1 x2) in
    let c = dxred env { c with x = c.x *. sgn mu1 } in
    holds "l2 casorati 2 is reduced" (c.ix = 0);
    rel_t "l2 casorati 2" 1.0e-9 c.x 1.0;
    let y1 = mul r.(i + 1) q.(i) in
    let y2 = mul r.(i) q.(i + 1) in
    let d = ref (sub y1 y2) in
    for j = 1 to (2 * mu1) - 1 do
      d :=
        dxadj env
          { !d with x = !d.x /. (nu +. dmu1 +. 1.0 -. float_of_int j) }
    done;
    if 2 * mu1 <= 1 then d := { !d with x = (nu +. 1.0) *. !d.x };
    let d = dxred env !d in
    holds "l2 casorati 1 is reduced" (d.ix = 0);
    rel_t "l2 casorati 1" 1.0e-9 d.x 1.0
  done

let casorati_mu dnu1 mu1 mu2 theta =
  let p = dxlegf env dnu1 0 mu1 mu2 theta Pneg in
  let q = dxlegf env dnu1 0 mu1 mu2 theta Q in
  let r = dxlegf env dnu1 0 mu1 mu2 theta Ppos in
  let sx = sin theta in
  for i = 0 to mu2 - mu1 - 1 do
    let mu = mu1 + i in
    let dmu = float_of_int mu in
    let x1 = mul p.(i + 1) q.(i) in
    let x2 = mul p.(i) q.(i + 1) in
    let x1 = { x1 with x = (dmu +. dnu1 +. 1.0) *. (dmu -. dnu1) *. x1.x } in
    let c = sub x1 x2 in
    let c = dxadj env { c with x = sx *. c.x *. sgn mu } in
    let c = dxred env c in
    holds "l2 casorati 4 is reduced" (c.ix = 0);
    rel_t "l2 casorati 4" 1.0e-9 c.x 1.0;
    let y1 = mul r.(i + 1) q.(i) in
    let y2 = mul r.(i) q.(i + 1) in
    let d = ref (sub y1 y2) in
    d := { !d with x = !d.x *. sx };
    for j = 1 to 2 * mu do
      d :=
        dxadj env
          { !d with x = !d.x /. (dnu1 +. dmu +. 1.0 -. float_of_int j) }
    done;
    let d = dxred env !d in
    holds "l2 casorati 3 is reduced" (d.ix = 0);
    rel_t "l2 casorati 3" 1.0e-9 d.x 1.0
  done

let l2_other_env () =
  let env2 = dxset 0 0 (1.0 /. min_float) 0 in
  Array.iter
    (fun v ->
      let b = dxred env2 (dxadj env2 { x = v; ix = 0 }) in
      holds "l2 round trip under an explicit dzero" (b.x = v && b.ix = 0))
    [| 1.0; -1.0; 0.5; -3.5; 123.456; 1.0e-300; 1.0e300; 0.1 |];
  holds "l2 an explicit dzero narrows the reduced band"
    ((dxred env2 (dxadj env2 { x = min_float; ix = 0 })).ix <> 0);
  let a = dxcon env { x = 1.0; ix = 5000 } in
  let b = dxcon env2 { x = 1.0; ix = 5000 } in
  holds "l2 dxcon agrees across environments" (a.x = b.x && a.ix = b.ix);
  let p = dxlegf env 100.0 0 10 10 0.9 Pnorm in
  let q = dxlegf env2 100.0 0 10 10 0.9 Pnorm in
  rel_t "l2 dxlegf agrees across environments" 1.0e-13 q.(0).x p.(0).x;
  let env3 = dxset 2 24 0.0 31 in
  let r = dxlegf env3 100.0 0 10 10 0.9 Pnorm in
  rel_t "l2 dxlegf with fewer series terms" 1.0e-5 r.(0).x p.(0).x

let l2_casorati () =
  casorati_nu 2.4 5 2 1.0;
  casorati_mu 2.4 2 7 1.0;
  casorati_nu 10.3 4 3 0.6;
  casorati_mu 10.3 3 8 0.6

(* ---- Level 3, the extended range ---- *)

let deg01 = 0.1 *. 4.0 *. atan 1.0 /. 180.0

let l3_casorati () =
  casorati_nu 2000.4 5 2000 deg01;
  casorati_mu 2000.4 1995 2000 deg01

let l3_nu_wise_against_mu_wise () =
  let dnu1 = 2000.4 and mu = 2000 in
  let a = dxlegf env dnu1 5 mu mu deg01 in
  let b = dxlegf env dnu1 0 (mu - 5) mu deg01 in
  List.iter
    (fun kind ->
      let u = dxcon env (a kind).(0) and v = dxcon env (b kind).(5) in
      rel_t "l3 nu-wise agrees with mu-wise" 1.0e-9 u.x v.x;
      holds "l3 nu-wise index agrees with mu-wise" (u.ix = v.ix))
    [ Pneg; Q; Ppos ]

let l3_extended_indices () =
  let p = dxlegf env 2000.4 5 2000 2000 deg01 Pneg in
  holds "l3 the index is used" (Array.for_all (fun a -> a.ix <> 0) p);
  let q = dxlegf env 2000.4 5 2000 2000 deg01 Q in
  holds "l3 Q needs the index the other way"
    (Array.for_all (fun a -> a.ix > 0) q);
  Array.iteri
    (fun i a ->
      let c = dxcon env a in
      holds "l3 dxcon normalises a huge value"
        (Float.abs c.x >= 0.1 && Float.abs c.x < 1.0);
      near_t "l3 dxcon keeps the value" 1.0e-6
        (log10 (Float.abs c.x) +. float_of_int c.ix)
        (log10 (Float.abs a.x) +. (float_of_int a.ix *. log10 2.0));
      ignore i)
    p

let l3_nrmp_against_legf () =
  List.iter
    (fun (n, m) ->
      let p, _ = dxnrmp env n m m deg01 By_theta in
      let q = dxlegf env (float_of_int n) 0 m m deg01 Pnorm in
      let u = dxcon env p.(0) and v = dxcon env q.(0) in
      rel_t "l3 dxnrmp agrees with dxlegf" 1.0e-11 u.x v.x;
      holds "l3 dxnrmp index agrees with dxlegf" (u.ix = v.ix))
    [ (100, 10); (100, 100); (500, 400); (1000, 1000) ]

(* ---- Domains ---- *)

let l1_domains () =
  let th = 1.0 in
  raises "dxlegf nudiff < 0 raises" (fun () ->
      ignore (dxlegf env 1.0 (-1) 0 0 th Pneg));
  raises "dxlegf dnu1 < -0.5 raises" (fun () ->
      ignore (dxlegf env (-0.6) 0 0 0 th Pneg));
  raises "dxlegf mu2 < mu1 raises" (fun () ->
      ignore (dxlegf env 1.0 0 3 2 th Pneg));
  raises "dxlegf mu1 < 0 raises" (fun () ->
      ignore (dxlegf env 1.0 0 (-1) 0 th Pneg));
  raises "dxlegf theta zero raises" (fun () ->
      ignore (dxlegf env 1.0 0 0 0 0.0 Pneg));
  raises "dxlegf theta above pi/2 raises" (fun () ->
      ignore (dxlegf env 1.0 0 0 0 2.0 Pneg));
  raises "dxlegf both ranges raises" (fun () ->
      ignore (dxlegf env 1.0 2 0 2 th Pneg));
  raises "dxlegf normalised P of fractional degree raises" (fun () ->
      ignore (dxlegf env 2.5 0 0 0 th Pnorm));
  raises "dxnrmp nu < 0 raises" (fun () ->
      ignore (dxnrmp env (-1) 0 0 0.5 By_x));
  raises "dxnrmp mu1 < 0 raises" (fun () ->
      ignore (dxnrmp env 2 (-1) 0 0.5 By_x));
  raises "dxnrmp mu1 > mu2 raises" (fun () ->
      ignore (dxnrmp env 2 3 2 0.5 By_x));
  raises "dxnrmp x out of range raises" (fun () ->
      ignore (dxnrmp env 2 0 0 1.5 By_x));
  raises "dxnrmp theta out of range raises" (fun () ->
      ignore (dxnrmp env 2 0 0 4.0 By_theta));
  raises "dxset bad radix raises" (fun () -> ignore (dxset 3 0 0.0 0));
  raises "dxset too few bits raises" (fun () -> ignore (dxset 0 0 0.0 14));
  raises "dxset too many bits raises" (fun () ->
      ignore (dxset 0 0 0.0 Sys.int_size));
  raises "dxset bad nradpl raises" (fun () -> ignore (dxset 0 200 0.0 0));
  raises "dxadj index overflow raises" (fun () ->
      ignore (dxadj env { x = 1.0; ix = max_int / 2 }));
  raises "dxcon index overflow raises" (fun () ->
      ignore (dxcon env { x = 1.0; ix = max_int / 2 }))

let l1_degenerate_vectors () =
  let th = 1.0 in
  let p = dxlegf env 3.0 0 5 5 th Ppos in
  holds "l1 P of positive order vanishes above the degree" (p.(0).x = 0.0);
  let p = dxlegf env 3.0 0 5 5 th Pnorm in
  holds "l1 normalised P vanishes above the degree" (p.(0).x = 0.0);
  let p = dxlegf env 3.0 0 2 5 th Pnorm in
  holds "l1 normalised P vanishes above the degree over mu"
    (p.(0).x <> 0.0 && p.(1).x <> 0.0 && p.(2).x = 0.0 && p.(3).x = 0.0)

let () =
  l1_roundtrip ();
  l1_zero ();
  l1_add_matches_float ();
  l1_extended ();
  l1_dxcon ();
  l1_legendre_closed ();
  l1_legendre_at_zero ();
  l1_nrmp_ends ();
  l1_degenerate_vectors ();
  l1_domains ();
  l2_arith ();
  l2_degree_recurrence ();
  l2_nrmp_modes ();
  l2_nrmp_against_legf ();
  l2_other_env ();
  l2_casorati ();
  l3_casorati ();
  l3_nu_wise_against_mu_wise ();
  l3_extended_indices ();
  l3_nrmp_against_legf ();
  if !bad = 0 then print_string "xlegf: all checks passed\n" else exit 1
