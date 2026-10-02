(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Level 1, level 2 and level 3 checks for the gamma and beta family, the
   polygamma sequence, Kummer's U and the Wigner 3j and 6j coefficients. *)

open Reeve.Gamma
open Reeve.Gamma.Raw

module G = Reeve.Gamma

let pi = 3.14159265358979323846
let euler_gamma = 0.5772156649015328606
let sq2pil = 0.91893853320467274178032973640562
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
let exact what a b = if not (a = b) then fail what

let rel_t what t a b =
  let e =
    if Float.abs b > 1.0e-14 then Float.abs (a -. b) /. Float.abs b
    else Float.abs (a -. b)
  in
  if not (e < t) then fail what

let holds what c = if not c then fail what

let raises what f =
  match f () with
  | _ -> fail what
  | exception Invalid_argument _ -> ()

(* ---- Exact integer arithmetic for the reference values ---- *)

let ibinom n k =
  let row = Array.make (n + 1) 0 in
  row.(0) <- 1;
  for i = 1 to n do
    for j = i downto 1 do
      row.(j) <- row.(j) + row.(j - 1)
    done
  done;
  row.(k)

let fact n =
  let p = ref 1.0 in
  for i = 1 to n do
    p := !p *. float_of_int i
  done;
  !p

(* The Racah closed form for a 3j symbol with all three projections zero. *)
let threej000 l1 l2 l3 =
  let s = l1 + l2 + l3 in
  if l1 + l2 < l3 || l2 + l3 < l1 || l3 + l1 < l2 then 0.0
  else if s land 1 = 1 then 0.0
  else begin
    let g = s / 2 in
    let num = fact (s - (2 * l1)) *. fact (s - (2 * l2)) *. fact (s - (2 * l3)) in
    let pre =
      fact g /. (fact (g - l1) *. fact (g - l2) *. fact (g - l3))
    in
    (if g land 1 = 0 then 1.0 else -1.0) *. sqrt (num /. fact (s + 1)) *. pre
  end

(* E1 from its own convergent series, for the complementary incomplete gamma
   at a = 0. *)
let e1 x =
  let s = ref 0.0 and term = ref 1.0 in
  let go = ref true and k = ref 1 in
  while !go && !k <= 200 do
    term := -. !term *. x;
    let t = !term /. (float_of_int !k *. fact !k) in
    s := !s +. t;
    if Float.abs t < 1.0e-18 *. Float.abs !s then go := false else incr k
  done;
  -.euler_gamma -. log x -. !s

(* ---- Level 1, regression ---- *)

let l1_factorial () =
  exact "l1 dfac 5" (dfac 5) 120.0;
  exact "l1 dfac 0" (dfac 0) 1.0

let l1_limits () =
  let xmin, xmax = dgamlm () in
  holds "l1 dgamlm straddles zero" (xmin < 0.0 && xmax > 0.0);
  near_t "l1 dgamlm xmax solves its equation" 1.0e-6
    (((xmax -. 0.5) *. log xmax) -. xmax +. 0.9189)
    (log max_float);
  near_t "l1 dgamlm xmin solves its equation" 1.0e-6
    ((((-.xmin) +. 0.5) *. log (-.xmin)) -. (-.xmin) -. 0.2258)
    (-.log min_float)

let l1_reciprocal () =
  for n = 0 to 20 do
    rel_t
      (Printf.sprintf "l1 dgamr %d" (n + 1))
      tol
      (dgamr (float_of_int (n + 1)))
      (1.0 /. fact n)
  done;
  exact "l1 dgamr 0" (dgamr 0.0) 0.0;
  exact "l1 dgamr -1" (dgamr (-1.0)) 0.0;
  exact "l1 dgamr -7" (dgamr (-7.0)) 0.0

let l1_incomplete () =
  rel_t "l1 dgami 1 2" tol (dgami 1.0 2.0) (1.0 -. exp (-2.0));
  rel_t "l1 dgamic 1 2" tol (dgamic 1.0 2.0) (exp (-2.0));
  rel_t "l1 dgamit 1 2" tol (dgamit 1.0 2.0) ((1.0 -. exp (-2.0)) /. 2.0);
  exact "l1 dgami a 0" (dgami 2.5 0.0) 0.0;
  rel_t "l1 dbetai 1 1" tol (dbetai 0.5 1.0 1.0) 0.5;
  exact "l1 dbetai at 0" (dbetai 0.0 2.0 3.0) 0.0;
  exact "l1 dbetai at 1" (dbetai 1.0 2.0 3.0) 1.0

(* ---- Level 2, mathematical identities ---- *)

let l2_log_gamma () =
  let lg x = fst (dlgams x) in
  let x = 3.5 in
  near "l2 dlgams recurrence" (lg (x +. 1.0)) (log x +. lg x);
  let x = 0.3 in
  loose "l2 dlgams reflection"
    (lg x +. lg (1.0 -. x))
    (log pi -. log (Float.abs (sin (pi *. x))));
  let x = 1.5 in
  loose "l2 dlgams duplication"
    (lg (2.0 *. x))
    (((2.0 *. x) -. 1.0) *. log 2.0 +. lg x +. lg (x +. 0.5)
    -. (0.5 *. log pi));
  holds "l2 dlgams sign below zero"
    (snd (dlgams (-0.5)) = -1.0
    && snd (dlgams (-1.5)) = 1.0
    && snd (dlgams (-2.5)) = -1.0
    && snd (dlgams 2.5) = 1.0);
  for n = 0 to 20 do
    let x = 10.0 +. (float_of_int n *. 7.5) in
    rel_t
      (Printf.sprintf "l2 d9lgmc at %.1f" x)
      1.0e-14 (lg x)
      ((((x -. 0.5) *. log x) -. x +. sq2pil) +. d9lgmc x)
  done

let l2_factorial () =
  for n = 0 to 20 do
    exact (Printf.sprintf "l2 dfac %d exact" n) (dfac n) (fact n)
  done;
  for n = 31 to 50 do
    rel_t
      (Printf.sprintf "l2 dfac %d recurrence" n)
      tol (dfac n)
      (float_of_int n *. dfac (n - 1))
  done;
  for n = 0 to 20 do
    for k = 0 to n do
      exact
        (Printf.sprintf "l2 dbinom %d %d exact" n k)
        (dbinom n k)
        (float_of_int (ibinom n k))
    done
  done;
  rel_t "l2 dbinom 60 30 symmetry" tol (dbinom 60 30)
    (dbinom 60 29 *. 31.0 /. 30.0)

let l2_poch () =
  exact "l2 dpoch a 0" (dpoch 2.5 0.0) 1.0;
  rel_t "l2 dpoch a 1" tol (dpoch 2.5 1.0) 2.5;
  rel_t "l2 dpoch 2.5 3" tol (dpoch 2.5 3.0) (2.5 *. 3.5 *. 4.5);
  for n = 0 to 20 do
    rel_t
      (Printf.sprintf "l2 dpoch 1 %d" n)
      tol
      (dpoch 1.0 (float_of_int n))
      (fact n)
  done;
  rel_t "l2 dpoch -5 3" tol (dpoch (-5.0) 3.0) (-60.0);
  rel_t "l2 dpoch -30 3" tol (dpoch (-30.0) 3.0) (-24360.0);
  exact "l2 dpoch -5 0.5" (dpoch (-5.0) 0.5) 0.0;
  let lg x = fst (dlgams x) and sg x = snd (dlgams x) in
  let a = 4.25 and x = 2.75 in
  rel_t "l2 dpoch against the gamma ratio" tol (dpoch a x)
    (sg (a +. x) *. sg a *. exp (lg (a +. x) -. lg a));
  let a = 30.5 and x = 21.25 in
  rel_t "l2 dpoch large against the gamma ratio" tol (dpoch a x)
    (sg (a +. x) *. sg a *. exp (lg (a +. x) -. lg a));
  let a = -7.25 and x = 0.5 in
  rel_t "l2 dpoch negative against the gamma ratio" loose_tol (dpoch a x)
    (sg (a +. x) *. sg a *. exp (lg (a +. x) -. lg a))

let l2_poch1 () =
  rel_t "l2 dpoch1 at zero is psi" tol (dpoch1 3.0 0.0) (1.5 -. euler_gamma);
  rel_t "l2 dpoch1 big x against dpoch" tol (dpoch1 3.0 2.0)
    ((dpoch 3.0 2.0 -. 1.0) /. 2.0);
  near_t "l2 dpoch1 tends to psi" 1.0e-7 (dpoch1 3.0 1.0e-8)
    (dpoch1 3.0 0.0);
  let step a x =
    rel_t
      (Printf.sprintf "l2 dpoch1 recurrence at %.2f" a)
      loose_tol
      (dpoch1 (a +. 1.0) x)
      ((1.0 +. ((a +. x) *. dpoch1 a x)) /. a)
  in
  step 5.0 1.0e-3;
  step 12.5 1.0e-3;
  step (-3.3) 1.0e-3;
  step (-8.7) 1.0e-4

let l2_beta () =
  for m = 1 to 6 do
    for n = 1 to 6 do
      let a = float_of_int m and b = float_of_int n in
      rel_t
        (Printf.sprintf "l2 dbeta %d %d" m n)
        tol (dbeta a b)
        (fact (m - 1) *. fact (n - 1) /. fact (m + n - 1))
    done
  done;
  rel_t "l2 dlbeta symmetry" tol (dlbeta 2.5 7.5) (dlbeta 7.5 2.5);
  rel_t "l2 dbeta against dlbeta" tol (dbeta 3.25 4.75)
    (exp (dlbeta 3.25 4.75));
  let step a b =
    rel_t
      (Printf.sprintf "l2 dlbeta recurrence at %.1f %.1f" a b)
      tol
      (dlbeta (a +. 1.0) b)
      (dlbeta a b +. log a -. log (a +. b))
  in
  step 0.5 0.75;
  step 2.5 7.5;
  step 12.0 25.0;
  step 100.0 100.0;
  step 200.0 300.0

let l2_betai () =
  let q = 3.0 in
  for i = 1 to 9 do
    let x = float_of_int i /. 10.0 in
    rel_t
      (Printf.sprintf "l2 dbetai p=1 at %.1f" x)
      tol (dbetai x 1.0 q)
      (1.0 -. ((1.0 -. x) ** q));
    rel_t
      (Printf.sprintf "l2 dbetai q=1 at %.1f" x)
      tol (dbetai x q 1.0) (x ** q);
    rel_t
      (Printf.sprintf "l2 dbetai symmetry at %.1f" x)
      tol (dbetai x 2.5 4.25)
      (1.0 -. dbetai (1.0 -. x) 4.25 2.5)
  done;
  holds "l2 dbetai increases"
    (dbetai 0.2 3.0 4.0 < dbetai 0.5 3.0 4.0
    && dbetai 0.5 3.0 4.0 < dbetai 0.8 3.0 4.0);
  rel_t "l2 dbetai large parameters" loose_tol
    (dbetai 0.5 200.0 200.0) 0.5

let l2_incomplete () =
  let lg x = fst (dlgams x) in
  for i = 1 to 6 do
    let a = float_of_int i /. 2.0 in
    for j = 1 to 6 do
      let x = float_of_int j in
      rel_t
        (Printf.sprintf "l2 dgami plus dgamic at %.1f %.1f" a x)
        tol
        (dgami a x +. dgamic a x)
        (exp (lg a));
      rel_t
        (Printf.sprintf "l2 dgamit scaling at %.1f %.1f" a x)
        tol (dgami a x)
        (exp (lg a +. (a *. log x)) *. dgamit a x)
    done
  done;
  let step a x =
    rel_t
      (Printf.sprintf "l2 dgamic recurrence at %.2f %.2f" a x)
      loose_tol
      (dgamic (a +. 1.0) x)
      ((a *. dgamic a x) +. ((x ** a) *. exp (-.x)))
  in
  step 2.5 1.5;
  step 0.5 0.25;
  step (-1.0) 0.5;
  step (-2.0) 0.5;
  step (-3.0) 0.5;
  step (-2.5) 0.75;
  for i = 1 to 4 do
    let x = float_of_int i /. 2.0 in
    rel_t
      (Printf.sprintf "l2 dgamic 0 is E1 at %.1f" x)
      loose_tol (dgamic 0.0 x) (e1 x)
  done;
  rel_t "l2 d9gmic 0 is E1" loose_tol (d9gmic 0.0 0.5 (log 0.5)) (e1 0.5);
  rel_t "l2 d9gmit 1 0.5" tol
    (d9gmit 1.0 0.5 (fst (dlgams 2.0)) 1.0)
    (2.0 *. (1.0 -. exp (-0.5)));
  rel_t "l2 d9lgic 1 2" tol (d9lgic 1.0 2.0 (log 2.0)) (-2.0);
  rel_t "l2 d9lgic 1 5" tol (d9lgic 1.0 5.0 (log 5.0)) (-5.0);
  rel_t "l2 d9lgit 2 2" tol
    (d9lgit 2.0 2.0 (log 2.0))
    (log ((1.0 -. (3.0 *. exp (-2.0))) /. 4.0))

let l2_polygamma () =
  rel_t "l2 dpsixn 1" hist_tol (dpsixn 1) (-.euler_gamma);
  rel_t "l2 dpsixn 2" hist_tol (dpsixn 2) (1.0 -. euler_gamma);
  for n = 1 to 150 do
    rel_t
      (Printf.sprintf "l2 dpsixn recurrence at %d" n)
      tol
      (dpsixn (n + 1))
      (dpsixn n +. (1.0 /. float_of_int n))
  done;
  let one x n kode =
    let a = Array.make 1 0.0 in
    let nz = dpsifn x n kode 1 a in
    exact (Printf.sprintf "l2 dpsifn nz zero at %.2f" x) (float_of_int nz) 0.0;
    a.(0)
  in
  rel_t "l2 dpsifn is minus psi at 1" hist_tol (one 1.0 0 1) euler_gamma;
  rel_t "l2 dpsifn is minus psi at 2" hist_tol (one 2.0 0 1)
    (euler_gamma -. 1.0);
  rel_t "l2 dpsifn kode 2" tol (one 2.0 0 2) (one 2.0 0 1 +. log 2.0);
  for k = 0 to 5 do
    List.iter
      (fun x ->
        rel_t
          (Printf.sprintf "l2 dpsifn recurrence k=%d at %.2f" k x)
          loose_tol (one x k 1)
          (one (x +. 1.0) k 1 +. (x ** float_of_int (-k - 1))))
      [ 0.25; 0.5; 1.0; 2.5; 7.0; 30.0 ]
  done;
  let m = 6 in
  let a = Array.make m 0.0 in
  let nz = dpsifn 1.75 2 1 m a in
  exact "l2 dpsifn sequence nz" (float_of_int nz) 0.0;
  for j = 0 to m - 1 do
    rel_t
      (Printf.sprintf "l2 dpsifn sequence member %d" j)
      tol a.(j)
      (one 1.75 (2 + j) 1)
  done;
  (* the scaled derivative at x = 1 is zeta(k+1), which tends to one *)
  for k = 20 to 120 do
    near_t
      (Printf.sprintf "l2 dpsifn is zeta at k=%d" k)
      1.0e-6 (one 1.0 k 1) 1.0
  done;
  (* the sequence degrades from the top when the top member overflows *)
  let m = 2000 in
  let a = Array.make m 0.0 in
  let nz = dpsifn 2.0 0 1 m a in
  holds "l2 dpsifn reports underflow" (nz > 0 && nz < m);
  rel_t "l2 dpsifn survives degradation" hist_tol a.(0) (euler_gamma -. 1.0);
  holds "l2 dpsifn zeroes the tail" (a.(m - 1) = 0.0 && a.(m - nz) = 0.0);
  let a = Array.make 3 0.0 in
  let nz = dpsifn 1.0e-20 0 1 3 a in
  exact "l2 dpsifn tiny x nz" (float_of_int nz) 0.0;
  exact "l2 dpsifn tiny x first" a.(0) 1.0e20;
  rel_t "l2 dpsifn tiny x second" tol a.(1) 1.0e40;
  rel_t "l2 dpsifn tiny x third" tol a.(2) 1.0e60

let l2_chu () =
  for i = 1 to 5 do
    let a = float_of_int i /. 2.0 in
    let b = a +. 1.0 in
    for j = 3 to 8 do
      let x = float_of_int j in
      rel_t
        (Printf.sprintf "l2 dchu a a+1 at %.1f %.1f" a x)
        tol (dchu a b x) (x ** -.a);
      rel_t
        (Printf.sprintf "l2 d9chu a a+1 at %.1f %.1f" a x)
        tol (d9chu a b x) 1.0
    done
  done;
  let kummer a b z =
    rel_t
      (Printf.sprintf "l2 dchu Kummer at %.2f %.2f %.2f" a b z)
      loose_tol (dchu a b z)
      ((z ** (1.0 -. b)) *. dchu (1.0 +. a -. b) (2.0 -. b) z)
  in
  kummer 0.3 0.7 0.5;
  kummer 0.5 1.4 0.1;
  kummer 1.25 0.4 0.75;
  kummer 0.6 1.3 0.5;
  let recur a b z =
    near_t
      (Printf.sprintf "l2 dchu recurrence at %.2f %.2f %.2f" a b z)
      loose_tol
      (dchu (a -. 1.0) b z
      +. ((b -. (2.0 *. a) -. z) *. dchu a b z)
      +. (a *. (a -. b +. 1.0) *. dchu (a +. 1.0) b z))
      0.0
  in
  recur 2.5 1.75 4.0;
  recur 1.5 0.5 3.0;
  recur 3.0 1.25 6.0;
  let ode a b z =
    let h = 1.0e-4 in
    let u = dchu a b z in
    let up = (dchu a b (z +. h) -. dchu a b (z -. h)) /. (2.0 *. h) in
    let upp =
      (dchu a b (z +. h) -. (2.0 *. u) +. dchu a b (z -. h)) /. (h *. h)
    in
    near_t
      (Printf.sprintf "l2 dchu solves its ode at %.2f %.2f %.2f" a b z)
      1.0e-4
      ((z *. upp) +. ((b -. z) *. up) -. (a *. u))
      0.0
  in
  ode 0.5 1.3 2.0;
  ode 1.5 2.5 3.0;
  holds "l2 dchu decays" (dchu 0.75 1.5 8.0 < dchu 0.75 1.5 4.0)

(* ---- The Wigner coefficients ---- *)

let steps lmin lmax = int_of_float (lmax -. lmin +. 1.5)

let l2_3jj () =
  let t = Array.make 40 0.0 in
  for i2 = 0 to 5 do
    for i3 = 0 to 5 do
      let l2 = float_of_int i2 and l3 = float_of_int i3 in
      let l1min, l1max = drc3jj l2 l3 0.0 0.0 t 40 in
      let n = steps l1min l1max in
      for j = 0 to n - 1 do
        let l1 = int_of_float (l1min +. float_of_int j +. 0.5) in
        rel_t
          (Printf.sprintf "l2 drc3jj (%d %d %d; 0 0 0)" l1 i2 i3)
          loose_tol t.(j)
          (threej000 l1 i2 i3)
      done;
      let s = ref 0.0 in
      for j = 0 to n - 1 do
        let l1 = l1min +. float_of_int j in
        s := !s +. ((l1 +. l1 +. 1.0) *. t.(j) *. t.(j))
      done;
      rel_t
        (Printf.sprintf "l2 drc3jj normalised %d %d" i2 i3)
        tol !s 1.0
    done
  done;
  (* half integral arguments, through the same normalisation *)
  let u = Array.make 40 0.0 in
  List.iter
    (fun (l2, l3, m2, m3) ->
      let l1min, l1max = drc3jj l2 l3 m2 m3 u 40 in
      let n = steps l1min l1max in
      let s = ref 0.0 in
      for j = 0 to n - 1 do
        let l1 = l1min +. float_of_int j in
        s := !s +. ((l1 +. l1 +. 1.0) *. u.(j) *. u.(j))
      done;
      rel_t
        (Printf.sprintf "l2 drc3jj normalised %.1f %.1f %.1f %.1f" l2 l3 m2 m3)
        tol !s 1.0)
    [
      (1.5, 1.0, 0.5, 0.0);
      (2.5, 1.5, -0.5, 0.5);
      (3.5, 3.5, 1.5, -2.5);
      (10.0, 8.0, 3.0, -2.0);
      (20.0, 17.0, 5.0, 4.0);
      (0.5, 0.5, 0.5, -0.5);
    ];
  (* the two 3j routines must agree where their ranges overlap *)
  let v = Array.make 40 0.0 in
  List.iter
    (fun (l2, l3, m2, m3) ->
      let l1min, l1max = drc3jj l2 l3 m2 m3 t 40 in
      let n = steps l1min l1max in
      for j = 0 to n - 1 do
        let l1 = l1min +. float_of_int j in
        let m1 = -.m2 -. m3 in
        let m2min, m2max = drc3jm l1 l2 l3 m1 v 40 in
        if m2 >= m2min -. 0.25 && m2 <= m2max +. 0.25 then begin
          let k = int_of_float (m2 -. m2min +. 0.5) in
          rel_t
            (Printf.sprintf "l2 drc3jj against drc3jm at l1=%.1f" l1)
            loose_tol v.(k) t.(j)
        end
      done;
      ignore l1max)
    [
      (2.0, 3.0, 1.0, -1.0);
      (1.5, 2.5, 0.5, 0.5);
      (4.0, 4.0, 2.0, -3.0);
      (6.0, 5.0, -2.0, 1.0);
    ]

let l2_3jm () =
  let t = Array.make 60 0.0 in
  List.iter
    (fun (l1, l2, l3, m1) ->
      let m2min, m2max = drc3jm l1 l2 l3 m1 t 60 in
      let n = steps m2min m2max in
      let s = ref 0.0 in
      for j = 0 to n - 1 do
        s := !s +. (t.(j) *. t.(j))
      done;
      rel_t
        (Printf.sprintf "l2 drc3jm orthogonal %.1f %.1f %.1f %.1f" l1 l2 l3 m1)
        tol !s
        (1.0 /. (l1 +. l1 +. 1.0)))
    [
      (1.0, 1.0, 1.0, 0.0);
      (2.0, 3.0, 4.0, 1.0);
      (1.5, 2.0, 2.5, 0.5);
      (3.5, 3.5, 3.0, -1.5);
      (10.0, 12.0, 8.0, 2.0);
      (20.0, 15.0, 18.0, -3.0);
    ];
  (* the j = 0 column, where the symbol has a closed form *)
  let m2min, m2max = drc3jm 2.0 2.0 0.0 1.0 t 60 in
  exact "l2 drc3jm single value range" (m2max -. m2min) 0.0;
  rel_t "l2 drc3jm (2 2 0; 1 -1 0)" tol t.(0)
    ((if (2 - 1) land 1 = 0 then 1.0 else -1.0) /. sqrt 5.0)

let l2_6j () =
  let t = Array.make 60 0.0 in
  (* the 6j symbol with a zero argument, where it has a closed form *)
  List.iter
    (fun (l2, l3) ->
      let l1min, l1max = drc6j l2 l3 0.0 l3 l2 t 60 in
      let n = steps l1min l1max in
      for j = 0 to n - 1 do
        let l1 = l1min +. float_of_int j in
        let phase = l1 +. l2 +. l3 in
        let sign =
          if Float.rem (Float.abs phase) 2.0 < 0.5 then 1.0 else -1.0
        in
        rel_t
          (Printf.sprintf "l2 drc6j zero argument %.1f %.1f at %.1f" l2 l3 l1)
          loose_tol t.(j)
          (sign /. sqrt (((l2 +. l2 +. 1.0) *. (l3 +. l3 +. 1.0))))
      done)
    [ (1.0, 1.0); (2.0, 1.0); (2.0, 2.0); (3.0, 2.0); (4.0, 4.0); (6.0, 5.0) ];
  (* orthogonality, which is what the normalisation rests on *)
  let u = Array.make 60 0.0 in
  List.iter
    (fun (l2, l3, l4, l5, l6) ->
      let l1min, l1max = drc6j l2 l3 l4 l5 l6 u 60 in
      let n = steps l1min l1max in
      let s = ref 0.0 in
      for j = 0 to n - 1 do
        let l1 = l1min +. float_of_int j in
        s := !s +. ((l1 +. l1 +. 1.0) *. (l4 +. l4 +. 1.0) *. u.(j) *. u.(j))
      done;
      rel_t
        (Printf.sprintf "l2 drc6j normalised %.1f %.1f %.1f %.1f %.1f" l2 l3 l4
           l5 l6)
        tol !s 1.0)
    [
      (1.5, 1.5, 1.0, 1.5, 1.5);
      (2.0, 3.0, 2.0, 3.0, 2.0);
      (3.5, 2.5, 3.0, 2.5, 3.5);
      (5.0, 4.0, 3.0, 4.0, 5.0);
      (10.0, 9.0, 8.0, 9.0, 10.0);
      (0.5, 0.5, 1.0, 0.5, 0.5);
    ];
  (* the 6j symbol is invariant under permuting its columns *)
  let v = Array.make 60 0.0 in
  List.iter
    (fun (l2, l3, l4, l5, l6) ->
      let l1min, _ = drc6j l2 l3 l4 l5 l6 u 60 in
      let l1min', l1max' = drc6j l3 l2 l4 l6 l5 v 60 in
      exact
        (Printf.sprintf "l2 drc6j column swap range %.1f %.1f" l2 l3)
        l1min l1min';
      let n = steps l1min' l1max' in
      for j = 0 to n - 1 do
        rel_t
          (Printf.sprintf "l2 drc6j column swap %.1f %.1f member %d" l2 l3 j)
          loose_tol v.(j) u.(j)
      done;
      let l1min'', l1max'' = drc6j l5 l6 l4 l2 l3 v 60 in
      exact
        (Printf.sprintf "l2 drc6j row swap range %.1f %.1f" l2 l3)
        l1min l1min'';
      let n = steps l1min'' l1max'' in
      for j = 0 to n - 1 do
        rel_t
          (Printf.sprintf "l2 drc6j row swap %.1f %.1f member %d" l2 l3 j)
          loose_tol v.(j) u.(j)
      done)
    [
      (2.0, 3.0, 2.0, 3.0, 2.0);
      (3.5, 2.5, 3.0, 2.5, 3.5);
      (4.0, 4.0, 4.0, 4.0, 4.0);
      (6.0, 5.0, 4.0, 5.0, 6.0);
    ]

(* ---- Level 3, classical tabulated values ---- *)

let l3_log_gamma () =
  let lg x = fst (dlgams x) in
  rel_t "l3 dlgams 1.5" hist_tol (lg 1.5) (-0.12078223763524522);
  rel_t "l3 dlgams 2.5" hist_tol (lg 2.5) 0.28468287047291918;
  rel_t "l3 dlgams 5" tol (lg 5.0) 3.1780538303479458

let l3_psi () =
  rel_t "l3 dpsixn 1" hist_tol (dpsixn 1) (-0.5772156649015329);
  rel_t "l3 dpsixn 2" hist_tol (dpsixn 2) 0.4227843350984671;
  let a = Array.make 1 0.0 in
  let _ = dpsifn 1.0 0 1 1 a in
  rel_t "l3 dpsifn psi 1" hist_tol (-.a.(0)) (-0.5772156649015329);
  let _ = dpsifn 2.0 0 1 1 a in
  rel_t "l3 dpsifn psi 2" hist_tol (-.a.(0)) 0.4227843350984671

(* ---- Domains the Fortran refuses ---- *)

let l1_domains () =
  raises "dlgams 0 raises" (fun () -> dlgams 0.0);
  raises "dlgams -3 raises" (fun () -> dlgams (-3.0));
  raises "d9lgmc below 10 raises" (fun () -> d9lgmc 9.5);
  raises "dfac negative raises" (fun () -> dfac (-1));
  raises "dfac 171 raises" (fun () -> dfac 171);
  raises "dbinom negative raises" (fun () -> dbinom (-1) 0);
  raises "dbinom n < m raises" (fun () -> dbinom 3 5);
  raises "dpoch pole raises" (fun () -> dpoch 0.5 (-2.5));
  raises "dlbeta non-positive raises" (fun () -> dlbeta 0.0 1.0);
  raises "dbeta non-positive raises" (fun () -> dbeta 1.0 (-1.0));
  raises "dbetai x outside raises" (fun () -> dbetai 1.5 1.0 1.0);
  raises "dbetai p non-positive raises" (fun () -> dbetai 0.5 0.0 1.0);
  raises "dgami a non-positive raises" (fun () -> dgami 0.0 1.0);
  raises "dgami x negative raises" (fun () -> dgami 1.0 (-1.0));
  raises "dgamic x negative raises" (fun () -> dgamic 1.0 (-1.0));
  raises "dgamic 0 0 raises" (fun () -> dgamic 0.0 0.0);
  raises "dgamit x negative raises" (fun () -> dgamit 1.0 (-1.0));
  raises "d9gmic a positive raises" (fun () -> d9gmic 1.0 0.5 (log 0.5));
  raises "d9gmic x non-positive raises" (fun () -> d9gmic (-1.0) 0.0 0.0);
  raises "d9gmit x non-positive raises" (fun () -> d9gmit 1.0 0.0 0.0 1.0);
  raises "d9lgit a below x raises" (fun () -> d9lgit 1.0 2.0 0.0);
  raises "d9lgit x non-positive raises" (fun () -> d9lgit 2.0 0.0 0.0);
  raises "dpsixn 0 raises" (fun () -> dpsixn 0);
  let a = Array.make 2 0.0 in
  raises "dpsifn x non-positive raises" (fun () -> dpsifn 0.0 0 1 1 a);
  raises "dpsifn n negative raises" (fun () -> dpsifn 1.0 (-1) 1 1 a);
  raises "dpsifn kode raises" (fun () -> dpsifn 1.0 0 3 1 a);
  raises "dpsifn m raises" (fun () -> dpsifn 1.0 0 1 0 a);
  raises "dpsifn short array raises" (fun () -> dpsifn 1.0 0 1 3 a);
  raises "dchu x zero raises" (fun () -> dchu 1.0 2.0 0.0);
  raises "dchu x negative raises" (fun () -> dchu 1.0 2.0 (-1.0));
  raises "dchu 1+a-b near zero raises" (fun () -> dchu 1.0 2.0 0.25);
  let t = Array.make 4 0.0 in
  raises "drc3jj l2 below m2 raises" (fun () ->
      drc3jj 1.0 2.0 2.0 0.0 t 4);
  raises "drc3jj non integral raises" (fun () ->
      drc3jj 1.3 2.0 0.0 0.0 t 4);
  raises "drc3jj array too small raises" (fun () ->
      drc3jj 6.0 6.0 0.0 0.0 t 4);
  raises "drc3jm triangle raises" (fun () -> drc3jm 1.0 1.0 5.0 0.0 t 4);
  raises "drc3jm l1 below m1 raises" (fun () -> drc3jm 1.0 2.0 2.0 2.0 t 4);
  raises "drc6j triangle raises" (fun () -> drc6j 1.0 1.0 9.0 1.0 1.0 t 4);
  raises "drc6j array too small raises" (fun () ->
      drc6j 6.0 6.0 6.0 6.0 6.0 t 4)

(* ---- The OCaml surface against the Fortran shapes it wraps ---- *)

let surface_psifn () =
  let m = 6 in
  let ar = Array.make m 0.0 in
  let nzr = dpsifn 1.75 2 1 m ar in
  let a, nz = G.dpsifn 1.75 2 m in
  exact "surface dpsifn nz" (float_of_int nz) (float_of_int nzr);
  holds "surface dpsifn length" (Array.length a = m);
  for j = 0 to m - 1 do
    exact (Printf.sprintf "surface dpsifn member %d" j) a.(j) ar.(j)
  done;
  let sr = Array.make 1 0.0 in
  ignore (dpsifn 2.0 0 2 1 sr);
  let s, _ = G.dpsifn ~scaled:true 2.0 0 1 in
  exact "surface dpsifn scaled" s.(0) sr.(0);
  let u, _ = G.dpsifn ~scaled:false 2.0 0 1 in
  rel_t "surface dpsifn unscaled" tol (s.(0) -. u.(0)) (log 2.0);
  let big = 2000 in
  let br = Array.make big 0.0 in
  let nzr = dpsifn 2.0 0 1 big br in
  let b, nz = G.dpsifn 2.0 0 big in
  holds "surface dpsifn reports underflow" (nz = nzr && nz > 0 && nz < big);
  holds "surface dpsifn length under underflow" (Array.length b = big);
  exact "surface dpsifn tail is zero" b.(big - 1) 0.0;
  for j = 0 to big - 1 do
    exact (Printf.sprintf "surface dpsifn long member %d" j) b.(j) br.(j)
  done;
  raises "surface dpsifn x non-positive raises" (fun () -> G.dpsifn 0.0 0 1);
  raises "surface dpsifn n negative raises" (fun () -> G.dpsifn 1.0 (-1) 1);
  raises "surface dpsifn m raises" (fun () -> G.dpsifn 1.0 0 0)

let surface_wigner () =
  let w = Array.make 64 0.0 in
  List.iter
    (fun (l2, l3, m2, m3) ->
      let lo, hi = drc3jj l2 l3 m2 m3 w 64 in
      let c, clo, chi = G.drc3jj l2 l3 m2 m3 in
      exact "surface drc3jj min" clo lo;
      exact "surface drc3jj max" chi hi;
      holds "surface drc3jj length" (Array.length c = steps lo hi);
      Array.iteri
        (fun i v ->
          exact (Printf.sprintf "surface drc3jj member %d" i) v w.(i))
        c)
    [
      (6.0, 4.0, -2.0, 1.0);
      (2.0, 3.0, 1.0, -1.0);
      (3.5, 3.5, 1.5, -2.5);
      (20.0, 17.0, 5.0, 4.0);
      (0.5, 0.5, 0.5, -0.5);
      (0.0, 3.0, 0.0, -1.0);
    ];
  List.iter
    (fun (l1, l2, l3, m1) ->
      let lo, hi = drc3jm l1 l2 l3 m1 w 64 in
      let c, clo, chi = G.drc3jm l1 l2 l3 m1 in
      exact "surface drc3jm min" clo lo;
      exact "surface drc3jm max" chi hi;
      holds "surface drc3jm length" (Array.length c = steps lo hi);
      Array.iteri
        (fun i v ->
          exact (Printf.sprintf "surface drc3jm member %d" i) v w.(i))
        c)
    [
      (1.0, 1.0, 1.0, 0.0);
      (2.0, 3.0, 4.0, 1.0);
      (1.5, 2.0, 2.5, 0.5);
      (3.5, 3.5, 3.0, -1.5);
      (2.0, 0.0, 2.0, 1.0);
    ];
  List.iter
    (fun (l2, l3, l4, l5, l6) ->
      let lo, hi = drc6j l2 l3 l4 l5 l6 w 64 in
      let c, clo, chi = G.drc6j l2 l3 l4 l5 l6 in
      exact "surface drc6j min" clo lo;
      exact "surface drc6j max" chi hi;
      holds "surface drc6j length" (Array.length c = steps lo hi);
      Array.iteri
        (fun i v -> exact (Printf.sprintf "surface drc6j member %d" i) v w.(i))
        c)
    [
      (1.5, 1.5, 1.0, 1.5, 1.5);
      (2.0, 3.0, 2.0, 3.0, 2.0);
      (10.0, 9.0, 8.0, 9.0, 10.0);
      (0.5, 0.5, 1.0, 0.5, 0.5);
      (0.0, 1.0, 1.0, 1.0, 1.0);
    ];
  raises "surface drc3jj l2 below m2 raises" (fun () ->
      G.drc3jj 1.0 2.0 2.0 0.0);
  raises "surface drc3jj non integral raises" (fun () ->
      G.drc3jj 1.3 2.0 0.0 0.0);
  raises "surface drc3jm triangle raises" (fun () -> G.drc3jm 1.0 1.0 5.0 0.0);
  raises "surface drc3jm l1 below m1 raises" (fun () ->
      G.drc3jm 1.0 2.0 2.0 2.0);
  raises "surface drc6j triangle raises" (fun () ->
      G.drc6j 1.0 1.0 9.0 1.0 1.0)

let () =
  l1_factorial ();
  l1_limits ();
  l1_reciprocal ();
  l1_incomplete ();
  l1_domains ();
  l2_log_gamma ();
  l2_factorial ();
  l2_poch ();
  l2_poch1 ();
  l2_beta ();
  l2_betai ();
  l2_incomplete ();
  l2_polygamma ();
  l2_chu ();
  l2_3jj ();
  l2_3jm ();
  l2_6j ();
  l3_log_gamma ();
  l3_psi ();
  surface_psifn ();
  surface_wigner ();
  if !bad = 0 then print_string "gamma: all checks passed\n" else exit 1
