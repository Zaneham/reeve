(* Reeve, Copyright 2026 Zane Hambly.
   The gamma and beta family: the reciprocal and log gamma, factorials and
   binomials, Pochhammer's symbol, the complete and incomplete beta, the
   incomplete gamma trio, the polygamma sequence, Kummer's confluent
   hypergeometric function U, and the Wigner 3j and 6j coefficients, ported
   from the SLATEC double precision routines. *)

open Mach

let pi = 3.141592653589793238462643383279503
let sq2pil = 0.91893853320467274178032973640562
let euler = 0.57721566490153286060651209008240

let ipow x n =
  let rec go acc b e =
    if e = 0 then acc
    else go (if e land 1 = 1 then acc *. b else acc) (b *. b) (e lsr 1)
  in
  if n < 0 then 1.0 /. go 1.0 x (-n) else go 1.0 x n

let parity n = if n land 1 = 0 then 1.0 else -1.0

(* ---- Coefficient tables ---- *)

let algmcs =
  [|
    0.1666389480451863247205729650822e+0;
    -0.1384948176067563840732986059135e-4;
    0.9810825646924729426157171547487e-8;
    -0.1809129475572494194263306266719e-10;
    0.6221098041892605227126015543416e-13;
    -0.3399615005417721944303330599666e-15;
    0.2683181998482698748957538846666e-17;
    -0.2868042435334643284144622399999e-19;
    0.3962837061046434803679306666666e-21;
    -0.6831888753985766870111999999999e-23;
    0.1429227355942498147573333333333e-24;
    -0.3547598158101070547199999999999e-26;
    0.1025680058010470912000000000000e-27;
    -0.3401102254316748799999999999999e-29;
    0.1276642195630062933333333333333e-30;
  |]

let facn =
  [|
    0.100000000000000000000000000000000e+1;
    0.100000000000000000000000000000000e+1;
    0.200000000000000000000000000000000e+1;
    0.600000000000000000000000000000000e+1;
    0.240000000000000000000000000000000e+2;
    0.120000000000000000000000000000000e+3;
    0.720000000000000000000000000000000e+3;
    0.504000000000000000000000000000000e+4;
    0.403200000000000000000000000000000e+5;
    0.362880000000000000000000000000000e+6;
    0.362880000000000000000000000000000e+7;
    0.399168000000000000000000000000000e+8;
    0.479001600000000000000000000000000e+9;
    0.622702080000000000000000000000000e+10;
    0.871782912000000000000000000000000e+11;
    0.130767436800000000000000000000000e+13;
    0.209227898880000000000000000000000e+14;
    0.355687428096000000000000000000000e+15;
    0.640237370572800000000000000000000e+16;
    0.121645100408832000000000000000000e+18;
    0.243290200817664000000000000000000e+19;
    0.510909421717094400000000000000000e+20;
    0.112400072777760768000000000000000e+22;
    0.258520167388849766400000000000000e+23;
    0.620448401733239439360000000000000e+24;
    0.155112100433309859840000000000000e+26;
    0.403291461126605635584000000000000e+27;
    0.108888694504183521607680000000000e+29;
    0.304888344611713860501504000000000e+30;
    0.884176199373970195454361600000000e+31;
    0.265252859812191058636308480000000e+33;
  |]

let bern =
  [|
    0.833333333333333333333333333333333e-1;
    -0.138888888888888888888888888888888e-2;
    0.330687830687830687830687830687830e-4;
    -0.826719576719576719576719576719576e-6;
    0.208767569878680989792100903212014e-7;
    -0.528419013868749318484768220217955e-9;
    0.133825365306846788328269809751291e-10;
    -0.338968029632258286683019539124944e-12;
    0.858606205627784456413590545042562e-14;
    -0.217486869855806187304151642386591e-15;
    0.550900282836022951520265260890225e-17;
    -0.139544646858125233407076862640635e-18;
    0.353470703962946747169322997780379e-20;
    -0.895351742703754685040261131811274e-22;
    0.226795245233768306031095073886816e-23;
    -0.574472439520264523834847971943400e-24;
    0.145517247561486490186626486727132e-26;
    -0.368599494066531017818178247990866e-28;
    0.933673425709504467203255515278562e-30;
    -0.236502241570062993455963519636983e-31;
  |]

let psifn_b =
  [|
    1.00000000000000000e+00;
    -5.00000000000000000e-01;
    1.66666666666666667e-01;
    -3.33333333333333333e-02;
    2.38095238095238095e-02;
    -3.33333333333333333e-02;
    7.57575757575757576e-02;
    -2.53113553113553114e-01;
    1.16666666666666667e+00;
    -7.09215686274509804e+00;
    5.49711779448621554e+01;
    -5.29124242424242424e+02;
    6.19212318840579710e+03;
    -8.65802531135531136e+04;
    1.42551716666666667e+06;
    -2.72982310678160920e+07;
    6.01580873900642368e+08;
    -1.51163157670921569e+10;
    4.29614643061166667e+11;
    -1.37116552050883328e+13;
    4.88332318973593167e+14;
    -1.92965793419400681e+16;
  |]

(* ---- Series lengths and derived constants ---- *)

let nalgm = 6

let d9lgmc_xbig = 1.0 /. sqrt eps_2_dp

let d9lgmc_xmax =
  exp (Float.min (log (huge_dp /. 12.0)) (-.log (12.0 *. tiny_dp)))

let dpoch1_sqtbig = 1.0 /. sqrt (24.0 *. tiny_dp)
let dpoch1_alneps = log eps_2_dp
let dbinom_bilnmx = log huge_dp -. 0.0001
let dbinom_fintmx = 0.9 /. eps_2_dp
let dbeta_xmax = 171.61447887182297
let dbeta_alnsml = log tiny_dp
let dbetai_eps = eps_2_dp
let dbetai_alneps = log eps_2_dp
let dbetai_sml = tiny_dp
let dbetai_alnsml = log tiny_dp
let gmic_eps = 0.5 *. eps_2_dp
let gmic_bot = log tiny_dp
let gamic_alneps = -.log eps_2_dp

(* ---- Elementary pieces the family needs ---- *)

let dpsi = Specfun.dpsi

(* There's no Fortran behind this one. DGAMMA and DLNGAM were dropped for the
   2008 intrinsics, so it's DGAMLN and reflection, and gam is the exp of it,
   which costs about 1e-13 by x = 170. *)
let lngam x =
  if x > 0.0 then Specfun.dgamln x
  else begin
    if Float.trunc x = x then
      invalid_arg "lngam: x is zero or a negative integer";
    let sinpix = Float.abs (sin (pi *. x)) in
    if sinpix = 0.0 then invalid_arg "lngam: x is a negative integer";
    log pi -. log sinpix -. Specfun.dgamln (1.0 -. x)
  end

let gamsgn x =
  if x > 0.0 then 1.0
  else begin
    let i = int_of_float (Float.rem (-.Float.trunc x) 2.0 +. 0.1) in
    if i = 0 then -1.0 else 1.0
  end

let gam x = gamsgn x *. exp (lngam x)

(* ---- Gamma limits, reciprocal and signed log ---- *)

(* The bound is 0.01 too wide. dgamlm.inc:64 reads
   Xmax = Xmax ! (Label 200 removed) - 0.01_DP, so the nudge the original
   applied is inside the comment now. Line 54 went the same way. *)
let dgamlm () =
  let alnsml = log tiny_dp and alnbig = log huge_dp in
  let xmin = ref (-.alnsml) and done_ = ref false in
  let i = ref 1 in
  while (not !done_) && !i <= 10 do
    let xold = !xmin in
    let xln = log !xmin in
    xmin :=
      !xmin
      -. (!xmin *. ((((!xmin +. 0.5) *. xln) -. !xmin -. 0.2258) +. alnsml)
          /. ((!xmin *. xln) +. 0.5));
    if Float.abs (!xmin -. xold) < 0.005 then done_ := true else incr i
  done;
  if not !done_ then invalid_arg "dgamlm: unable to find xmin";
  let xmin = ref (-. !xmin) in
  let xmax = ref alnbig and done_ = ref false in
  let i = ref 1 in
  while (not !done_) && !i <= 10 do
    let xold = !xmax in
    let xln = log !xmax in
    xmax :=
      !xmax
      -. (!xmax *. ((((!xmax -. 0.5) *. xln) -. !xmax +. 0.9189) -. alnbig)
          /. ((!xmax *. xln) -. 0.5));
    if Float.abs (!xmax -. xold) < 0.005 then done_ := true else incr i
  done;
  if not !done_ then invalid_arg "dgamlm: unable to find xmax";
  xmin := Float.max !xmin (-. !xmax +. 1.0);
  (!xmin, !xmax)

let dlgams x =
  let dlgam = lngam x in
  if x > 0.0 then (dlgam, 1.0)
  else begin
    let i = int_of_float (Float.rem (-.Float.trunc x) 2.0 +. 0.1) in
    if i = 0 then (dlgam, -1.0) else (dlgam, 1.0)
  end

let dgamr x =
  if x <= 0.0 && Float.trunc x = x then 0.0
  else if Float.abs x > 10.0 then begin
    let alngx, sgngx = dlgams x in
    sgngx *. exp (-.alngx)
  end
  else 1.0 /. gam x

let d9lgmc x =
  if x < 10.0 then invalid_arg "d9lgmc: x must be >= 10"
  else if x < d9lgmc_xbig then begin
    let t = 10.0 /. x in
    Specfun.dcsevl ((2.0 *. t *. t) -. 1.0) algmcs nalgm /. x
  end
  else if x < d9lgmc_xmax then 1.0 /. (12.0 *. x)
  else 0.0

(* ---- Factorials, binomials and Pochhammer ---- *)

let dfac n =
  if n < 0 then invalid_arg "dfac: factorial of a negative integer undefined"
  else if n <= 30 then facn.(n)
  else if n <= 170 then begin
    let x = float_of_int n +. 1.0 in
    exp (((x -. 0.5) *. log x) -. x +. sq2pil +. d9lgmc x)
  end
  else invalid_arg "dfac: n so big that factorial n overflows"

let dbinom n m =
  if n < 0 || m < 0 then invalid_arg "dbinom: n or m < 0";
  if n < m then invalid_arg "dbinom: n < m";
  let k = min m (n - m) in
  let small =
    if k <= 20 && float_of_int k *. log (float_of_int (max n 1)) <= dbinom_bilnmx
    then begin
      let b = ref 1.0 in
      if k = 0 then Some 1.0
      else begin
        for i = 1 to k do
          let xn = float_of_int (n - i + 1) and xk = float_of_int i in
          b := !b *. (xn /. xk)
        done;
        Some (if !b < dbinom_fintmx then Float.trunc (!b +. 0.5) else !b)
      end
    end
    else None
  in
  match small with
  | Some v -> v
  | None ->
    if k < 9 then
      invalid_arg "dbinom: result overflows because n and/or m too big";
    let xn = float_of_int (n + 1) in
    let xk = float_of_int (k + 1) in
    let xnk = float_of_int (n - k + 1) in
    let corr = d9lgmc xn -. d9lgmc xk -. d9lgmc xnk in
    let r =
      (xk *. log (xnk /. xk))
      -. (xn *. Expint.dlnrel (-.(xk -. 1.0) /. xn))
      -. (0.5 *. log (xn *. xnk /. xk))
      +. 1.0 -. sq2pil +. corr
    in
    if r > dbinom_bilnmx then
      invalid_arg "dbinom: result overflows because n and/or m too big";
    let r = exp r in
    if r < dbinom_fintmx then Float.trunc (r +. 0.5) else r

let dpoch a x =
  let ax = a +. x in
  if ax <= 0.0 && Float.trunc ax = ax then begin
    if a > 0.0 || Float.trunc a <> a then
      invalid_arg "dpoch: a+x is a non-positive integer but a is not";
    if x = 0.0 then 1.0
    else begin
      let n = int_of_float x in
      if Float.min (a +. x) a < -20.0 then
        parity n
        *. exp
             (((a -. 0.5) *. Expint.dlnrel (x /. (a -. 1.0)))
             +. (x *. log (-.a +. 1.0 -. x))
             -. x +. d9lgmc (-.a +. 1.0)
             -. d9lgmc (-.a -. x +. 1.0))
      else begin
        let ia = int_of_float a in
        parity n *. dfac (-ia) /. dfac (-ia - n)
      end
    end
  end
  else if a <= 0.0 && Float.trunc a = a then 0.0
  else begin
    let n = int_of_float (Float.abs x) in
    let general =
      if float_of_int n <> x || n > 20 then begin
        let absax = Float.abs (a +. x) and absa = Float.abs a in
        if Float.max absax absa <= 20.0 then Some (gam (a +. x) *. dgamr a)
        else if Float.abs x > 0.5 *. absa then begin
          let alngax, sgngax = dlgams (a +. x) in
          let alnga, sgnga = dlgams a in
          Some (sgngax *. sgnga *. exp (alngax -. alnga))
        end
        else None
      end
      else begin
        let p = ref 1.0 in
        for i = 1 to n do
          p := !p *. (a +. float_of_int (i - 1))
        done;
        Some !p
      end
    in
    match general with
    | Some v -> v
    | None ->
      let b = if a < 0.0 then -.a -. x +. 1.0 else a in
      let r =
        exp
          (((b -. 0.5) *. Expint.dlnrel (x /. b))
          +. (x *. log (b +. x))
          -. x +. d9lgmc (b +. x) -. d9lgmc b)
      in
      if a < 0.0 && r <> 0.0 then
        r /. (cos (pi *. x) +. (Specfun.dcot (pi *. a) *. sin (pi *. x)))
      else r
  end

let dpoch1 a x =
  if x = 0.0 then dpsi a
  else begin
    let absx = Float.abs x and absa = Float.abs a in
    if absx > 0.1 *. absa then (dpoch a x -. 1.0) /. x
    else if absx *. log (Float.max absa 2.0) > 0.1 then (dpoch a x -. 1.0) /. x
    else begin
      let bp = if a < -0.5 then 1.0 -. a -. x else a in
      let incr = if bp < 10.0 then 11 - int_of_float bp else 0 in
      let b = bp +. float_of_int incr in
      let var = b +. (0.5 *. (x -. 1.0)) in
      let alnvar = log var in
      let q = x *. alnvar in
      let poly1 = ref 0.0 in
      if var < dpoch1_sqtbig then begin
        let vinv = 1.0 /. var in
        let var2 = vinv *. vinv in
        let rho = 0.5 *. (x +. 1.0) in
        let gbern = Array.make 21 0.0 in
        gbern.(0) <- 1.0;
        gbern.(1) <- -.rho /. 12.0;
        let term = ref var2 in
        poly1 := gbern.(1) *. !term;
        let nterms = int_of_float (-0.5 *. dpoch1_alneps /. alnvar) + 1 in
        if nterms > 20 then
          invalid_arg "dpoch1: nterms is too big, maybe eps is bad";
        if nterms >= 2 then
          for k = 2 to nterms do
            let gbk = ref 0.0 in
            for j = 1 to k do
              let ndx = k - j + 1 in
              gbk := !gbk +. (bern.(ndx - 1) *. gbern.(j - 1))
            done;
            gbern.(k) <- -.rho *. !gbk /. float_of_int k;
            let fk = float_of_int k in
            term :=
              !term *. ((2.0 *. fk) -. 2.0 -. x) *. ((2.0 *. fk) -. 1.0 -. x)
              *. var2;
            poly1 := !poly1 +. (gbern.(k) *. !term)
          done
      end;
      poly1 := (x -. 1.0) *. !poly1;
      let r = ref ((Expint.dexprl q *. (alnvar +. (q *. !poly1))) +. !poly1) in
      if incr <> 0 then
        for ii = 1 to incr do
          let i = incr - ii in
          let binv = 1.0 /. (bp +. float_of_int i) in
          r := (!r -. binv) /. (1.0 +. (x *. binv))
        done;
      if bp = a then !r
      else begin
        let sinpxx = sin (pi *. x) /. x in
        let sinpx2 = sin (0.5 *. pi *. x) in
        let trig =
          (sinpxx *. Specfun.dcot (pi *. b)) -. (2.0 *. sinpx2 *. (sinpx2 /. x))
        in
        trig +. ((1.0 +. (x *. trig)) *. !r)
      end
    end
  end

(* ---- Beta ---- *)

let dlbeta a b =
  let p = Float.min a b and q = Float.max a b in
  if p <= 0.0 then invalid_arg "dlbeta: both arguments must be > 0"
  else if p >= 10.0 then begin
    let corr = d9lgmc p +. d9lgmc q -. d9lgmc (p +. q) in
    (-0.5 *. log q) +. sq2pil +. corr
    +. ((p -. 0.5) *. log (p /. (p +. q)))
    +. (q *. Expint.dlnrel (-.p /. (p +. q)))
  end
  else if q < 10.0 then log (gam p *. (gam q /. gam (p +. q)))
  else begin
    let corr = d9lgmc q -. d9lgmc (p +. q) in
    lngam p +. corr +. p -. (p *. log (p +. q))
    +. ((q -. 0.5) *. Expint.dlnrel (-.p /. (p +. q)))
  end

let dbeta a b =
  if a <= 0.0 || b <= 0.0 then invalid_arg "dbeta: both arguments must be > 0";
  if a +. b < dbeta_xmax then gam a *. gam b /. gam (a +. b)
  else begin
    let r = dlbeta a b in
    if r < dbeta_alnsml then 0.0 else exp r
  end

let dbetai x pin qin =
  if x < 0.0 || x > 1.0 then invalid_arg "dbetai: x is not in the range (0,1)";
  if pin <= 0.0 || qin <= 0.0 then invalid_arg "dbetai: p and/or q is <= 0";
  let y = ref x and p = ref pin and q = ref qin in
  if (qin > pin || x >= 0.8) && x >= 0.2 then begin
    y := 1.0 -. !y;
    p := qin;
    q := pin
  end;
  let y = !y and p = !p and q = !q in
  if (p +. q) *. y /. (p +. 1.0) < dbetai_eps then begin
    let xb = (p *. log (Float.max y dbetai_sml)) -. log p -. dlbeta p q in
    let r = if xb > dbetai_alnsml && y <> 0.0 then exp xb else 0.0 in
    if y <> x || p <> pin then 1.0 -. r else r
  end
  else begin
    let ps = q -. Float.trunc q in
    let ps = if ps = 0.0 then 1.0 else ps in
    let xb = (p *. log y) -. dlbeta ps p -. log p in
    let acc = ref 0.0 in
    if xb >= dbetai_alnsml then begin
      acc := exp xb;
      let term = ref (!acc *. p) in
      if ps <> 1.0 then begin
        let n = int_of_float (Float.max (dbetai_alneps /. log y) 4.0) in
        for i = 1 to n do
          let xi = float_of_int i in
          term := !term *. (xi -. ps) *. y /. xi;
          acc := !acc +. (!term /. (p +. xi))
        done
      end
    end;
    if q > 1.0 then begin
      let xb = (p *. log y) +. (q *. log (1.0 -. y)) -. dlbeta p q -. log q in
      let ib = ref (int_of_float (Float.max (xb /. dbetai_alnsml) 0.0)) in
      let term = ref (exp (xb -. (float_of_int !ib *. dbetai_alnsml))) in
      let c = 1.0 /. (1.0 -. y) in
      let p1 = q *. c /. (p +. q -. 1.0) in
      let finsum = ref 0.0 in
      let n = int_of_float q in
      let n = if q = float_of_int n then n - 1 else n in
      let go = ref true in
      let i = ref 1 in
      while !go && !i <= n do
        if p1 <= 1.0 && !term /. dbetai_eps <= !finsum then go := false
        else begin
          let xi = float_of_int !i in
          term := (q -. xi +. 1.0) *. c *. !term /. (p +. q -. xi);
          if !term > 1.0 then begin
            ib := !ib - 1;
            term := !term *. dbetai_sml
          end;
          if !ib = 0 then finsum := !finsum +. !term;
          incr i
        end
      done;
      acc := !acc +. !finsum
    end;
    let r = if y <> x || p <> pin then 1.0 -. !acc else !acc in
    Float.max (Float.min r 1.0) 0.0
  end

(* ---- Incomplete gamma ---- *)

let d9lgic a x alx =
  let xpa = x +. 1.0 -. a in
  let xma = x -. 1.0 -. a in
  let r = ref 0.0 and p = ref 1.0 and s = ref 1.0 in
  let go = ref true and k = ref 1 in
  while !go && !k <= 300 do
    let fk = float_of_int !k in
    let t = fk *. (a -. fk) *. (1.0 +. !r) in
    r := -.t /. (((xma +. (2.0 *. fk)) *. (xpa +. (2.0 *. fk))) +. t);
    p := !r *. !p;
    s := !s +. !p;
    if Float.abs !p < gmic_eps *. !s then go := false else incr k
  done;
  if !go then
    invalid_arg "d9lgic: no convergence in 300 terms of the continued fraction";
  (a *. alx) -. x +. log (!s /. xpa)

let d9lgit a x algap1 =
  if x <= 0.0 || a < x then invalid_arg "d9lgit: x should be > 0 and <= a";
  let ax = a +. x in
  let a1x = ax +. 1.0 in
  let r = ref 0.0 and p = ref 1.0 and s = ref 1.0 in
  let go = ref true and k = ref 1 in
  while !go && !k <= 200 do
    let fk = float_of_int !k in
    let t = (a +. fk) *. x *. (1.0 +. !r) in
    r := t /. (((ax +. fk) *. (a1x +. fk)) -. t);
    p := !r *. !p;
    s := !s +. !p;
    if Float.abs !p < gmic_eps *. !s then go := false else incr k
  done;
  if !go then
    invalid_arg "d9lgit: no convergence in 200 terms of the continued fraction";
  let hstar = 1.0 -. (x *. !s /. a1x) in
  -.x -. algap1 -. log hstar

let d9gmic a x alx =
  if a > 0.0 then invalid_arg "d9gmic: a must be near a negative integer";
  if x <= 0.0 then invalid_arg "d9gmic: x must be > 0";
  let m = int_of_float (-.(a -. 0.5)) in
  let fm = float_of_int m in
  let te = ref 1.0 and t = ref 1.0 and s = ref 1.0 in
  let go = ref true and k = ref 1 in
  while !go && !k <= 200 do
    let fkp1 = float_of_int (!k + 1) in
    te := -.x *. !te /. (fm +. fkp1);
    t := !te /. fkp1;
    s := !s +. !t;
    if Float.abs !t < gmic_eps *. !s then go := false else incr k
  done;
  if !go then
    invalid_arg "d9gmic: no convergence in 200 terms of the continued fraction";
  let r = ref (-.alx -. euler +. (x *. !s /. (fm +. 1.0))) in
  if m = 0 then !r
  else if m = 1 then -. !r -. 1.0 +. (1.0 /. x)
  else begin
    let te = ref fm and t = ref 1.0 and s = ref 1.0 in
    let go = ref true and k = ref 1 in
    while !go && !k <= m - 1 do
      let fk = float_of_int !k in
      te := -.x *. !te /. fk;
      t := !te /. (fm -. fk);
      s := !s +. !t;
      if Float.abs !t < gmic_eps *. Float.abs !s then go := false else incr k
    done;
    for k = 1 to m do
      r := !r +. (1.0 /. float_of_int k)
    done;
    let sgng = if m land 1 = 1 then -1.0 else 1.0 in
    let alng = log !r -. lngam (fm +. 1.0) in
    let r = ref (if alng > gmic_bot then sgng *. exp alng else 0.0) in
    if !s <> 0.0 then
      r :=
        !r
        +. Float.copy_sign (exp ((-.fm *. alx) +. log (Float.abs !s /. fm))) !s;
    !r
  end

let d9gmit a x algap1 sgngam =
  if x <= 0.0 then invalid_arg "d9gmit: x should be > 0";
  let ma = if a < 0.0 then int_of_float (a -. 0.5) else int_of_float (a +. 0.5) in
  let aeps = a -. float_of_int ma in
  let ae = if a < -0.5 then aeps else a in
  let t = ref 1.0 and te = ref ae and s = ref 1.0 in
  let go = ref true and k = ref 1 in
  while !go && !k <= 200 do
    let fk = float_of_int !k in
    te := -.x *. !te /. fk;
    t := !te /. (ae +. fk);
    s := !s +. !t;
    if Float.abs !t < gmic_eps *. Float.abs !s then go := false else incr k
  done;
  if !go then
    invalid_arg "d9gmit: no convergence in 200 terms of the Taylor series";
  if a >= -0.5 then exp (-.algap1 +. log !s)
  else begin
    let algs = ref (-.lngam (1.0 +. aeps) +. log !s) in
    let s = ref 1.0 in
    let m = -ma - 1 in
    if m <> 0 then begin
      let t = ref 1.0 in
      let go = ref true and k = ref 1 in
      while !go && !k <= m do
        t := x *. !t /. (aeps -. float_of_int (m + 1 - !k));
        s := !s +. !t;
        if Float.abs !t < gmic_eps *. Float.abs !s then go := false else incr k
      done
    end;
    algs := (-.float_of_int ma *. log x) +. !algs;
    if !s = 0.0 || aeps = 0.0 then exp !algs
    else begin
      let sgng2 = sgngam *. Float.copy_sign 1.0 !s in
      let alg2 = -.x -. algap1 +. log (Float.abs !s) in
      if alg2 > gmic_bot then (sgng2 *. exp alg2) +. exp !algs else 0.0
    end
  end


let dgamit a x =
  if x < 0.0 then invalid_arg "dgamit: x is negative";
  let alx = if x <> 0.0 then log x else 0.0 in
  let sga = if a <> 0.0 then Float.copy_sign 1.0 a else 1.0 in
  let ainta = Float.trunc (a +. (0.5 *. sga)) in
  let aeps = a -. ainta in
  if x = 0.0 then
    if ainta > 0.0 || aeps <> 0.0 then dgamr (a +. 1.0) else 0.0
  else if x <= 1.0 then begin
    let algap1, sgngam =
      if a >= -0.5 || aeps <> 0.0 then dlgams (a +. 1.0) else (0.0, 1.0)
    in
    d9gmit a x algap1 sgngam
  end
  else if a < x then begin
    let alng = d9lgic a x alx in
    let h = ref 1.0 and out = ref None in
    if aeps <> 0.0 || ainta > 0.0 then begin
      let algap1, sgngam = dlgams (a +. 1.0) in
      let t = log (Float.abs a) +. alng -. algap1 in
      if t > gamic_alneps then
        out := Some (-.sga *. sgngam *. exp (t -. (a *. alx)))
      else if t > -.gamic_alneps then h := 1.0 -. (sga *. sgngam *. exp t)
    end;
    match !out with
    | Some v -> v
    | None ->
      let t = (-.a *. alx) +. log (Float.abs !h) in
      Float.copy_sign (exp t) !h
  end
  else exp (d9lgit a x (lngam (a +. 1.0)))

let dgamic a x =
  if x < 0.0 then invalid_arg "dgamic: x is negative";
  if x > 0.0 then begin
    let alx = log x in
    let sga = if a <> 0.0 then Float.copy_sign 1.0 a else 1.0 in
    let ainta = Float.trunc (a +. (0.5 *. sga)) in
    let aeps = a -. ainta in
    let early = ref None and izero = ref false in
    let algap1 = ref 0.0 and sgngam = ref 1.0 in
    let sgngs = ref 1.0 and alngs = ref 0.0 in
    if x >= 1.0 then begin
      if a < x then early := Some (exp (d9lgic a x alx))
      else begin
        sgngam := 1.0;
        algap1 := lngam (a +. 1.0);
        sgngs := 1.0;
        alngs := d9lgit a x !algap1
      end
    end
    else begin
      let quick =
        if a <= 0.5 && Float.abs aeps <= 0.001 then begin
          let e =
            if -.ainta > 1.0 then
              2.0 *. (-.ainta +. 2.0) /. ((ainta *. ainta) -. 1.0)
            else 2.0
          in
          let e = e -. (alx *. (x ** -0.001)) in
          if e *. Float.abs aeps <= gmic_eps then Some (d9gmic a x alx) else None
        end
        else None
      in
      match quick with
      | Some v -> early := Some v
      | None ->
        let g, s = dlgams (a +. 1.0) in
        algap1 := g;
        sgngam := s;
        let gstar = d9gmit a x !algap1 !sgngam in
        if gstar = 0.0 then izero := true
        else begin
          alngs := log (Float.abs gstar);
          sgngs := Float.copy_sign 1.0 gstar
        end
    end;
    match !early with
    | Some v -> v
    | None ->
      let h = ref 1.0 and out = ref None in
      if not !izero then begin
        let t = (a *. alx) +. !alngs in
        if t > gamic_alneps then begin
          let sgng = -. !sgngs *. sga *. !sgngam in
          let t = t +. !algap1 -. log (Float.abs a) in
          out := Some (sgng *. exp t)
        end
        else if t > -.gamic_alneps then h := 1.0 -. (!sgngs *. exp t)
      end;
      (match !out with
       | Some v -> v
       | None ->
         let sgng = Float.copy_sign 1.0 !h *. sga *. !sgngam in
         let t = log (Float.abs !h) +. !algap1 -. log (Float.abs a) in
         sgng *. exp t)
  end
  else if a <= 0.0 then
    invalid_arg "dgamic: x = 0 and a <= 0 so dgamic is undefined"
  else exp (lngam (a +. 1.0) -. log a)

let dgami a x =
  if a <= 0.0 then invalid_arg "dgami: a must be > 0";
  if x < 0.0 then invalid_arg "dgami: x must be >= 0"
  else if x = 0.0 then 0.0
  else begin
    let factor = exp (lngam a +. (a *. log x)) in
    factor *. dgamit a x
  end

(* ---- Polygamma ---- *)

let dpsixn = Expint.dpsixn

let psifn_nmax = 100

(* This is off the F77, because dpsifn.inc types everything REAL(SP) and then
   fills its tables with _DP literals. *)
let dpsifn x n kode m ans =
  if x <= 0.0 then invalid_arg "dpsifn: x <= 0";
  if n < 0 then invalid_arg "dpsifn: n < 0";
  if kode < 1 || kode > 2 then invalid_arg "dpsifn: kode is neither 1 nor 2";
  if m < 1 then invalid_arg "dpsifn: m < 1";
  if Array.length ans < m then invalid_arg "dpsifn: ans is shorter than m";
  let nz = ref 0 and mm = ref m in
  let elim =
    2.302
    *. ((float_of_int (min (-min_exp_dp) max_exp_dp) *. log10_radix_dp) -. 3.0)
  in
  let wdtol = Float.max (eps_dp *. 0.5) 0.5e-18 in
  let xln = log x in
  let trm = Array.make 22 0.0 and trmr = Array.make psifn_nmax 0.0 in
  let s = ref 0.0 and xdmy = ref x and xdmln = ref xln and xinc = ref 0.0 in
  let nx = ref 0 in
  let n0rec = ref false and tail = ref false in
  let returned = ref false and brk = ref false in
  while not (!brk || !returned) do
    let nn = n + !mm - 1 in
    let fn = ref (float_of_int nn) in
    let fnp = ref (!fn +. 1.0) in
    let t = !fnp *. xln in
    if Float.abs t <= elim then
      if x < wdtol then begin
        ans.(0) <- ipow x (-n - 1);
        if !mm <> 1 then
          for i = 2 to !mm do
            ans.(i - 1) <- ans.(i - 2) /. x
          done;
        if n = 0 && kode = 2 then ans.(0) <- ans.(0) +. xln;
        returned := true
      end
      else begin
        let rln = Float.min (log10_radix_dp *. float_of_int digits_dp) 18.06 in
        let fln = Float.max rln 3.0 -. 3.0 in
        let yint = 3.50 +. (0.40 *. fln) in
        let slope = 0.21 +. (fln *. ((0.0006038 *. fln) +. 0.008677)) in
        let xm = yint +. (slope *. !fn) in
        let xmin = float_of_int (int_of_float xm + 1) in
        let series = ref false in
        let fln = ref fln in
        if n <> 0 then begin
          let xm = (-2.302 *. rln) -. Float.min 0.0 xln in
          let arg = Float.min 0.0 (xm /. float_of_int n) in
          let eps = exp arg in
          let xm = if Float.abs arg < 1.0e-3 then -.arg else 1.0 -. eps in
          fln := x *. xm /. eps;
          if xmin -. x > 7.0 && !fln < 15.0 then series := true
        end;
        if !series then begin
          let nn = int_of_float !fln + 1 in
          let np = n + 1 in
          let t = ref (exp (-.(float_of_int (n + 1) *. xln))) in
          s := !t;
          let den = ref x in
          for i = 1 to nn do
            den := !den +. 1.0;
            trm.(i - 1) <- ipow !den (-np);
            s := !s +. trm.(i - 1)
          done;
          ans.(0) <- !s;
          if n = 0 && kode = 2 then ans.(0) <- !s +. xln;
          if !mm <> 1 then begin
            let tol = wdtol /. 5.0 in
            for j = 2 to !mm do
              t := !t /. x;
              s := !t;
              let tols = !t *. tol in
              let den = ref x in
              let go = ref true and i = ref 1 in
              while !go && !i <= nn do
                den := !den +. 1.0;
                trm.(!i - 1) <- trm.(!i - 1) /. !den;
                s := !s +. trm.(!i - 1);
                if trm.(!i - 1) < tols then go := false else incr i
              done;
              ans.(j - 1) <- !s
            done
          end;
          returned := true
        end
        else begin
          xdmy := x;
          xdmln := xln;
          xinc := 0.0;
          if x < xmin then begin
            nx := int_of_float x;
            xinc := xmin -. float_of_int !nx;
            xdmy := x +. !xinc;
            xdmln := log !xdmy
          end;
          let t = !fn *. !xdmln in
          let t1 = !xdmln +. !xdmln in
          let t2 = t +. !xdmln in
          let tk = Float.max (Float.abs t) (Float.max (Float.abs t1) (Float.abs t2)) in
          if tk > elim then begin
            incr nz;
            ans.(!mm - 1) <- 0.0;
            decr mm;
            if !mm = 0 then returned := true
          end
          else begin
            let tss = ref (exp (-.t)) in
            let tt = 0.5 /. !xdmy in
            let tst = wdtol *. tt in
            let t1 = ref tt in
            if nn <> 0 then t1 := tt +. (1.0 /. !fn);
            let rxsq = 1.0 /. (!xdmy *. !xdmy) in
            let ta = 0.5 *. rxsq in
            let t = ref (!fnp *. ta) in
            s := !t *. psifn_b.(2);
            if Float.abs !s >= tst then begin
              let tk = ref 2.0 in
              let go = ref true and k = ref 4 in
              while !go && !k <= 22 do
                t :=
                  !t
                  *. ((!tk +. !fn +. 1.0) /. (!tk +. 1.0))
                  *. ((!tk +. !fn) /. (!tk +. 2.0))
                  *. rxsq;
                trm.(!k - 1) <- !t *. psifn_b.(!k - 1);
                if Float.abs trm.(!k - 1) < tst then go := false
                else begin
                  s := !s +. trm.(!k - 1);
                  tk := !tk +. 2.0;
                  incr k
                end
              done
            end;
            s := (!s +. !t1) *. !tss;
            let bad = ref false in
            if !xinc <> 0.0 then begin
              nx := int_of_float !xinc;
              let np = nn + 1 in
              if !nx > psifn_nmax then begin
                nz := 0;
                invalid_arg "dpsifn: n too large for the backward recurrence"
              end
              else if nn = 0 then begin
                n0rec := true;
                brk := true;
                bad := true
              end
              else begin
                let xm = ref (!xinc -. 1.0) in
                let fx = ref (x +. !xm) in
                for i = 1 to !nx do
                  trmr.(i - 1) <- ipow !fx (-np);
                  s := !s +. trmr.(i - 1);
                  xm := !xm -. 1.0;
                  fx := x +. !xm
                done
              end
            end;
            if not !bad then begin
              ans.(!mm - 1) <- !s;
              if !fn = 0.0 then begin
                tail := true;
                brk := true
              end
              else if !mm = 1 then returned := true
              else begin
                let j = ref 2 in
                let stop = ref false in
                while (not !stop) && !j <= !mm do
                  fnp := !fn;
                  fn := !fn -. 1.0;
                  tss := !tss *. !xdmy;
                  t1 := tt;
                  if !fn <> 0.0 then t1 := tt +. (1.0 /. !fn);
                  let t = ref (!fnp *. ta) in
                  s := !t *. psifn_b.(2);
                  if Float.abs !s >= tst then begin
                    let tk = ref (3.0 +. !fnp) in
                    let go = ref true and k = ref 4 in
                    while !go && !k <= 22 do
                      trm.(!k - 1) <- trm.(!k - 1) *. !fnp /. !tk;
                      if Float.abs trm.(!k - 1) < tst then go := false
                      else begin
                        s := !s +. trm.(!k - 1);
                        tk := !tk +. 2.0;
                        incr k
                      end
                    done
                  end;
                  s := (!s +. !t1) *. !tss;
                  let skip = ref false in
                  if !xinc <> 0.0 then
                    if !fn = 0.0 then begin
                      n0rec := true;
                      brk := true;
                      stop := true;
                      skip := true
                    end
                    else begin
                      let xm = ref (!xinc -. 1.0) in
                      let fx = ref (x +. !xm) in
                      for i = 1 to !nx do
                        trmr.(i - 1) <- trmr.(i - 1) *. !fx;
                        s := !s +. trmr.(i - 1);
                        xm := !xm -. 1.0;
                        fx := x +. !xm
                      done
                    end;
                  if not !skip then begin
                    ans.(!mm - !j) <- !s;
                    if !fn = 0.0 then begin
                      tail := true;
                      brk := true;
                      stop := true
                    end
                    else incr j
                  end
                done;
                if not !stop then returned := true
              end
            end
          end
        end
      end
    else if t <= 0.0 then begin
      nz := 0;
      invalid_arg "dpsifn: overflow, x too small or n+m-1 too large or both"
    end
    else begin
      incr nz;
      ans.(!mm - 1) <- 0.0;
      decr mm;
      if !mm = 0 then returned := true
    end
  done;
  if !n0rec || !tail then begin
    (* Grouped as the Fortran groups it, (X+NX)-I, which cancels where doing it
       in ints wouldn't. *)
    if !n0rec then
      for i = 1 to !nx do
        s := !s +. (1.0 /. (x +. float_of_int !nx -. float_of_int i))
      done;
    if kode = 2 then begin
      if !xdmy <> x then ans.(0) <- !s -. log (!xdmy /. x)
    end
    else ans.(0) <- !s -. !xdmln
  end;
  !nz

(* ---- Kummer's confluent hypergeometric function U ---- *)

let d9chu a b z =
  let eps = 4.0 *. eps_dp in
  let aa = Array.make 4 0.0 and bb = Array.make 4 0.0 in
  let bp = 1.0 +. a -. b in
  let ab = a *. bp in
  let ct2 = ref (2.0 *. (z -. ab)) in
  let sab = a +. bp in
  bb.(0) <- 1.0;
  aa.(0) <- 1.0;
  let ct3 = ref (sab +. 1.0 +. ab) in
  bb.(1) <- 1.0 +. (2.0 *. z /. !ct3);
  aa.(1) <- 1.0 +. (!ct2 /. !ct3);
  let anbn = ref (!ct3 +. sab +. 3.0) in
  let ct1 = ref (1.0 +. (2.0 *. z /. !anbn)) in
  bb.(2) <- 1.0 +. (6.0 *. !ct1 *. z /. !ct3);
  aa.(2) <- 1.0 +. (6.0 *. ab /. !anbn) +. (3.0 *. !ct1 *. !ct2 /. !ct3);
  let go = ref true and i = ref 4 in
  while !go && !i <= 300 do
    let x2i1 = float_of_int ((2 * !i) - 3) in
    ct1 := x2i1 /. (x2i1 -. 2.0);
    anbn := !anbn +. x2i1 +. sab;
    ct2 := (x2i1 -. 1.0) /. !anbn;
    let c2 = (x2i1 *. !ct2) -. 1.0 in
    let d1z = x2i1 *. 2.0 *. z /. !anbn in
    ct3 := sab *. !ct2;
    let g1 = d1z +. (!ct1 *. (c2 +. !ct3)) in
    let g2 = d1z -. c2 in
    let g3 = !ct1 *. (1.0 -. !ct3 -. (2.0 *. !ct2)) in
    bb.(3) <- (g1 *. bb.(2)) +. (g2 *. bb.(1)) +. (g3 *. bb.(0));
    aa.(3) <- (g1 *. aa.(2)) +. (g2 *. aa.(1)) +. (g3 *. aa.(0));
    if Float.abs ((aa.(3) *. bb.(0)) -. (aa.(0) *. bb.(3)))
       < eps *. Float.abs (bb.(3) *. bb.(0))
    then go := false
    else begin
      for j = 1 to 3 do
        aa.(j - 1) <- aa.(j);
        bb.(j - 1) <- bb.(j)
      done;
      incr i
    end
  done;
  if !go then invalid_arg "d9chu: no convergence in 300 terms";
  aa.(3) /. bb.(3)

let dchu a b x =
  if x = 0.0 then invalid_arg "dchu: x is zero so dchu is infinite";
  if x < 0.0 then invalid_arg "dchu: x is negative, use cchu";
  let eps = eps_2_dp in
  if Float.max (Float.abs a) 1.0 *. Float.max (Float.abs (1.0 +. a -. b)) 1.0
     < 0.99 *. Float.abs x
  then (x ** -.a) *. d9chu a b x
  else begin
    if Float.abs (1.0 +. a -. b) < sqrt eps then
      invalid_arg "dchu: the algorithm is bad when 1+a-b is near zero for small x";
    let aintb = if b >= 0.0 then Float.trunc (b +. 0.5) else Float.trunc (b -. 0.5) in
    let beps = b -. aintb in
    let n = int_of_float aintb in
    let alnx = log x in
    let xtoeps = exp (-.beps *. alnx) in
    let summ = ref 0.0 in
    if n >= 1 then begin
      let m = n - 2 in
      if m >= 0 then begin
        let t = ref 1.0 in
        summ := 1.0;
        if m <> 0 then
          for i = 1 to m do
            let xi = float_of_int i in
            t := !t *. (a -. b +. xi) *. x /. ((1.0 -. b +. xi) *. xi);
            summ := !summ +. !t
          done;
        summ := gam (b -. 1.0) *. dgamr a *. ipow x (1 - n) *. xtoeps *. !summ
      end
    end
    else begin
      summ := 1.0;
      if n <> 0 then begin
        let t = ref 1.0 in
        let m = -n in
        for i = 1 to m do
          let xi1 = float_of_int (i - 1) in
          t := !t *. (a +. xi1) *. x /. ((b +. xi1) *. (xi1 +. 1.0));
          summ := !summ +. !t
        done
      end;
      summ := dpoch (1.0 +. a -. b) (-.a) *. !summ
    end;
    let istrt = if n < 1 then 1 - n else 0 in
    let xi = float_of_int istrt in
    let factor = parity n *. dgamr (1.0 +. a -. b) *. ipow x istrt in
    let factor =
      if beps <> 0.0 then factor *. beps *. pi /. sin (beps *. pi) else factor
    in
    let pochai = dpoch a xi in
    let gamri1 = dgamr (xi +. 1.0) in
    let gamrni = dgamr (aintb +. xi) in
    let b0 = ref (factor *. dpoch a (xi -. beps) *. gamrni *. dgamr (xi +. 1.0 -. beps)) in
    if Float.abs (xtoeps -. 1.0) <= 0.5 then begin
      let pch1ai = dpoch1 (a +. xi) (-.beps) in
      let pch1i = dpoch1 (xi +. 1.0 -. beps) beps in
      let c0 =
        ref
          (factor *. pochai *. gamrni *. gamri1
          *. (-.dpoch1 (b +. xi) (-.beps) +. pch1ai -. pch1i
             +. (beps *. pch1ai *. pch1i)))
      in
      let xeps1 = alnx *. Expint.dexprl (-.beps *. alnx) in
      let acc = ref (!summ +. !c0 +. (xeps1 *. !b0)) in
      let xn = float_of_int n in
      let go = ref true and i = ref 1 in
      while !go && !i <= 1000 do
        let xi = float_of_int (istrt + !i) in
        let xi1 = float_of_int (istrt + !i - 1) in
        b0 := (a +. xi1 -. beps) *. !b0 *. x /. ((xn +. xi1) *. (xi -. beps));
        c0 :=
          ((a +. xi1) *. !c0 *. x /. ((b +. xi1) *. xi))
          -. (((((a -. 1.0) *. (xn +. (2.0 *. xi) -. 1.0)) +. (xi *. (xi -. beps)))
              *. !b0)
             /. (xi *. (b +. xi1) *. (a +. xi1 -. beps)));
        let t = !c0 +. (xeps1 *. !b0) in
        acc := !acc +. t;
        if Float.abs t < eps *. Float.abs !acc then go := false else incr i
      done;
      if !go then
        invalid_arg "dchu: no convergence in 1000 terms of the ascending series";
      !acc
    end
    else begin
      let a0 = ref (factor *. pochai *. dgamr (b +. xi) *. gamri1 /. beps) in
      b0 := xtoeps *. !b0 /. beps;
      let acc = ref (!summ +. !a0 -. !b0) in
      let go = ref true and i = ref 1 in
      while !go && !i <= 1000 do
        let xi = float_of_int (istrt + !i) in
        let xi1 = float_of_int (istrt + !i - 1) in
        a0 := (a +. xi1) *. !a0 *. x /. ((b +. xi1) *. xi);
        b0 := (a +. xi1 -. beps) *. !b0 *. x /. ((aintb +. xi1) *. (xi -. beps));
        let t = !a0 -. !b0 in
        acc := !acc +. t;
        if Float.abs t < eps *. Float.abs !acc then go := false else incr i
      done;
      if !go then
        invalid_arg "dchu: no convergence in 1000 terms of the ascending series";
      !acc
    end
  end

(* ---- Wigner 3j and 6j coefficients ---- *)

let wigner_eps = 0.01

let drc3jj l2 l3 m2 m3 thrcof ndim =
  if ndim < 1 then invalid_arg "drc3jj: ndim < 1";
  if Array.length thrcof < ndim then
    invalid_arg "drc3jj: thrcof is shorter than ndim";
  let eps = wigner_eps in
  let hugee = sqrt (huge_dp /. 20.0) in
  let srhuge = sqrt hugee in
  let tinyy = 1.0 /. hugee in
  let srtiny = 1.0 /. srhuge in
  let m1 = -.m2 -. m3 in
  if l2 -. Float.abs m2 +. eps < 0.0 || l3 -. Float.abs m3 +. eps < 0.0 then
    invalid_arg "drc3jj: l2-abs m2 or l3-abs m3 less than zero";
  if Float.rem (l2 +. Float.abs m2 +. eps) 1.0 >= eps +. eps
     || Float.rem (l3 +. Float.abs m3 +. eps) 1.0 >= eps +. eps
  then invalid_arg "drc3jj: l2+abs m2 or l3+abs m3 not integer";
  let l1min = Float.max (Float.abs (l2 -. l3)) (Float.abs m1) in
  let l1max = l2 +. l3 in
  if Float.rem (l1max -. l1min +. eps) 1.0 >= eps +. eps then
    invalid_arg "drc3jj: l1max-l1min not integer";
  if not (l1min < l1max -. eps) then begin
    if not (l1min < l1max +. eps) then invalid_arg "drc3jj: l1min > l1max";
    thrcof.(0) <-
      parity (int_of_float (Float.abs (l2 +. m2 -. l3 +. m3) +. eps))
      /. sqrt (l1min +. l2 +. l3 +. 1.0);
    (l1min, l1max)
  end
  else begin
    let nfin = int_of_float (l1max -. l1min +. 1.0 +. eps) in
    if ndim < nfin then
      invalid_arg "drc3jj: result array for the 3j coefficients too small";
    let l1 = ref l1min in
    let newfac = ref 0.0 and oldfac = ref 0.0 in
    let c1 = ref 0.0 and c1old = ref 0.0 and denom = ref 0.0 in
    thrcof.(0) <- srtiny;
    let sum1 = ref ((!l1 +. !l1 +. 1.0) *. tinyy) in
    let sum2 = ref 0.0 and sumfor = ref 0.0 and sumbac = ref 0.0 in
    let sumuni = ref 0.0 in
    let lstep = ref 1 in
    let brk = ref false in
    while not !brk do
      incr lstep;
      l1 := !l1 +. 1.0;
      oldfac := !newfac;
      let a1 =
        (!l1 +. l2 +. l3 +. 1.0) *. (!l1 -. l2 +. l3) *. (!l1 +. l2 -. l3)
        *. (-. !l1 +. l2 +. l3 +. 1.0)
      in
      let a2 = (!l1 +. m1) *. (!l1 -. m1) in
      newfac := sqrt (a1 *. a2);
      if !l1 < 1.0 +. eps then
        c1 := -.(!l1 +. !l1 -. 1.0) *. !l1 *. (m3 -. m2) /. !newfac
      else begin
        let dv =
          (-.l2 *. (l2 +. 1.0) *. m1)
          +. (l3 *. (l3 +. 1.0) *. m1)
          +. (!l1 *. (!l1 -. 1.0) *. (m3 -. m2))
        in
        denom := (!l1 -. 1.0) *. !newfac;
        if !lstep > 2 then c1old := Float.abs !c1;
        c1 := -.(!l1 +. !l1 -. 1.0) *. dv /. !denom
      end;
      if !lstep <= 2 then begin
        let x = srtiny *. !c1 in
        thrcof.(1) <- x;
        sum1 := !sum1 +. (tinyy *. (!l1 +. !l1 +. 1.0) *. !c1 *. !c1);
        if !lstep = nfin then begin
          sumuni := !sum1;
          brk := true
        end
      end
      else begin
        let c2 = -. !l1 *. !oldfac /. !denom in
        let x = ref ((!c1 *. thrcof.(!lstep - 2)) +. (c2 *. thrcof.(!lstep - 3))) in
        thrcof.(!lstep - 1) <- !x;
        sumfor := !sum1;
        sum1 := !sum1 +. ((!l1 +. !l1 +. 1.0) *. !x *. !x);
        let cycle = ref false in
        if !lstep <> nfin then begin
          if Float.abs !x >= srhuge then begin
            for i = 1 to !lstep do
              if Float.abs thrcof.(i - 1) < srtiny then thrcof.(i - 1) <- 0.0;
              thrcof.(i - 1) <- thrcof.(i - 1) /. srhuge
            done;
            sum1 := !sum1 /. hugee;
            sumfor := !sumfor /. hugee;
            x := !x /. srhuge
          end;
          if !c1old > Float.abs !c1 then cycle := true
        end;
        if not !cycle then begin
          let x1 = !x in
          let x2 = thrcof.(!lstep - 2) in
          let x3 = thrcof.(!lstep - 3) in
          let nstep2 = nfin - !lstep + 3 in
          let nfinp2 = nfin + 2 in
          let nfinp3 = nfin + 3 in
          l1 := l1max;
          thrcof.(nfin - 1) <- srtiny;
          sum2 := tinyy *. (!l1 +. !l1 +. 1.0);
          l1 := !l1 +. 2.0;
          let lstep = ref 1 in
          let y = ref 0.0 in
          let bbrk = ref false in
          while not !bbrk do
            incr lstep;
            l1 := !l1 -. 1.0;
            oldfac := !newfac;
            let a1s =
              (!l1 +. l2 +. l3) *. (!l1 -. l2 +. l3 -. 1.0)
              *. (!l1 +. l2 -. l3 -. 1.0)
              *. (-. !l1 +. l2 +. l3 +. 2.0)
            in
            let a2s = (!l1 +. m1 -. 1.0) *. (!l1 -. m1 -. 1.0) in
            newfac := sqrt (a1s *. a2s);
            let dv =
              (-.l2 *. (l2 +. 1.0) *. m1)
              +. (l3 *. (l3 +. 1.0) *. m1)
              +. (!l1 *. (!l1 -. 1.0) *. (m3 -. m2))
            in
            denom := !l1 *. !newfac;
            c1 := -.(!l1 +. !l1 -. 1.0) *. dv /. !denom;
            if !lstep <= 2 then begin
              y := srtiny *. !c1;
              thrcof.(nfin - 2) <- !y;
              sumbac := !sum2;
              sum2 := !sum2 +. (tinyy *. (!l1 +. !l1 -. 3.0) *. !c1 *. !c1)
            end
            else begin
              let c2 = -.(!l1 -. 1.0) *. !oldfac /. !denom in
              y :=
                (!c1 *. thrcof.(nfinp2 - !lstep - 1))
                +. (c2 *. thrcof.(nfinp3 - !lstep - 1));
              if !lstep = nstep2 then begin
                let y3 = !y in
                let y2 = thrcof.(nfinp2 - !lstep - 1) in
                let y1 = thrcof.(nfinp3 - !lstep - 1) in
                let ratio =
                  ((x1 *. y1) +. (x2 *. y2) +. (x3 *. y3))
                  /. ((x1 *. x1) +. (x2 *. x2) +. (x3 *. x3))
                in
                let nlim = nfin - nstep2 + 1 in
                if Float.abs ratio < 1.0 then begin
                  let nlim = nlim + 1 in
                  let ratio = 1.0 /. ratio in
                  for nn = nlim to nfin do
                    thrcof.(nn - 1) <- ratio *. thrcof.(nn - 1)
                  done;
                  sumuni := !sumfor +. (ratio *. ratio *. !sumbac)
                end
                else begin
                  for nn = 1 to nlim do
                    thrcof.(nn - 1) <- ratio *. thrcof.(nn - 1)
                  done;
                  sumuni := (ratio *. ratio *. !sumfor) +. !sumbac
                end;
                bbrk := true
              end
              else begin
                thrcof.(nfin + 1 - !lstep - 1) <- !y;
                sumbac := !sum2;
                sum2 := !sum2 +. ((!l1 +. !l1 -. 3.0) *. !y *. !y);
                if Float.abs !y >= srhuge then begin
                  for i = 1 to !lstep do
                    let ix = nfin - i + 1 in
                    if Float.abs thrcof.(ix - 1) < srtiny then thrcof.(ix - 1) <- 0.0;
                    thrcof.(ix - 1) <- thrcof.(ix - 1) /. srhuge
                  done;
                  sum2 := !sum2 /. hugee;
                  sumbac := !sumbac /. hugee
                end
              end
            end
          done;
          brk := true
        end
      end
    done;
    let cnorm = ref (1.0 /. sqrt !sumuni) in
    let sign1 = Float.copy_sign 1.0 thrcof.(nfin - 1) in
    let sign2 = parity (int_of_float (Float.abs (l2 +. m2 -. l3 +. m3) +. eps)) in
    if sign1 *. sign2 <= 0.0 then cnorm := -. !cnorm;
    if Float.abs !cnorm < 1.0 then begin
      let thresh = tinyy /. Float.abs !cnorm in
      for nn = 1 to nfin do
        if Float.abs thrcof.(nn - 1) < thresh then thrcof.(nn - 1) <- 0.0;
        thrcof.(nn - 1) <- !cnorm *. thrcof.(nn - 1)
      done
    end
    else
      for nn = 1 to nfin do
        thrcof.(nn - 1) <- !cnorm *. thrcof.(nn - 1)
      done;
    (l1min, l1max)
  end

let drc3jm l1 l2 l3 m1 thrcof ndim =
  if ndim < 1 then invalid_arg "drc3jm: ndim < 1";
  if Array.length thrcof < ndim then
    invalid_arg "drc3jm: thrcof is shorter than ndim";
  let eps = wigner_eps in
  let hugee = sqrt (huge_dp /. 20.0) in
  let srhuge = sqrt hugee in
  let tinyy = 1.0 /. hugee in
  let srtiny = 1.0 /. srhuge in
  if l1 -. Float.abs m1 +. eps < 0.0
     || Float.rem (l1 +. Float.abs m1 +. eps) 1.0 >= eps +. eps
  then invalid_arg "drc3jm: l1-abs m1 less than zero or l1+abs m1 not integer";
  if l1 +. l2 -. l3 < -.eps || l1 -. l2 +. l3 < -.eps || -.l1 +. l2 +. l3 < -.eps
  then invalid_arg "drc3jm: l1, l2, l3 do not satisfy the triangular condition";
  if Float.rem (l1 +. l2 +. l3 +. eps) 1.0 >= eps +. eps then
    invalid_arg "drc3jm: l1+l2+l3 not integer";
  let m2min = Float.max (-.l2) (-.l3 -. m1) in
  let m2max = Float.min l2 (l3 -. m1) in
  if Float.rem (m2max -. m2min +. eps) 1.0 >= eps +. eps then
    invalid_arg "drc3jm: m2max-m2min not integer";
  if not (m2min < m2max -. eps) then begin
    if not (m2min < m2max +. eps) then invalid_arg "drc3jm: m2min > m2max";
    thrcof.(0) <-
      parity (int_of_float (Float.abs (l2 -. l3 -. m1) +. eps))
      /. sqrt (l1 +. l2 +. l3 +. 1.0);
    (m2min, m2max)
  end
  else begin
    let nfin = int_of_float (m2max -. m2min +. 1.0 +. eps) in
    if ndim < nfin then
      invalid_arg "drc3jm: result array for the 3j coefficients too small";
    let m2 = ref m2min in
    thrcof.(0) <- srtiny;
    let newfac = ref 0.0 and oldfac = ref 0.0 in
    let c1 = ref 0.0 and c1old = ref 0.0 in
    let sum1 = ref tinyy and sum2 = ref 0.0 in
    let sumfor = ref 0.0 and sumbac = ref 0.0 and sumuni = ref 0.0 in
    let lstep = ref 1 in
    let brk = ref false in
    while not !brk do
      incr lstep;
      m2 := !m2 +. 1.0;
      let m3 = -.m1 -. !m2 in
      oldfac := !newfac;
      let a1 = (l2 -. !m2 +. 1.0) *. (l2 +. !m2) *. (l3 +. m3 +. 1.0) *. (l3 -. m3) in
      newfac := sqrt a1;
      let dv =
        ((l1 +. l2 +. l3 +. 1.0) *. (l2 +. l3 -. l1))
        -. ((l2 -. !m2 +. 1.0) *. (l3 +. m3 +. 1.0))
        -. ((l2 +. !m2 -. 1.0) *. (l3 -. m3 -. 1.0))
      in
      if !lstep > 2 then c1old := Float.abs !c1;
      c1 := -.dv /. !newfac;
      if !lstep <= 2 then begin
        let x = srtiny *. !c1 in
        thrcof.(1) <- x;
        sum1 := !sum1 +. (tinyy *. !c1 *. !c1);
        if !lstep = nfin then begin
          sumuni := !sum1;
          brk := true
        end
      end
      else begin
        let c2 = -. !oldfac /. !newfac in
        let x = ref ((!c1 *. thrcof.(!lstep - 2)) +. (c2 *. thrcof.(!lstep - 3))) in
        thrcof.(!lstep - 1) <- !x;
        sumfor := !sum1;
        sum1 := !sum1 +. (!x *. !x);
        let cycle = ref false in
        if !lstep <> nfin then begin
          if Float.abs !x >= srhuge then begin
            for i = 1 to !lstep do
              if Float.abs thrcof.(i - 1) < srtiny then thrcof.(i - 1) <- 0.0;
              thrcof.(i - 1) <- thrcof.(i - 1) /. srhuge
            done;
            sum1 := !sum1 /. hugee;
            sumfor := !sumfor /. hugee;
            x := !x /. srhuge
          end;
          if !c1old > Float.abs !c1 then cycle := true
        end;
        if not !cycle then begin
          let nstep2 = nfin - !lstep + 3 in
          let x1 = !x in
          let x2 = thrcof.(!lstep - 2) in
          let x3 = thrcof.(!lstep - 3) in
          let nfinp2 = nfin + 2 in
          let nfinp3 = nfin + 3 in
          thrcof.(nfin - 1) <- srtiny;
          sum2 := tinyy;
          m2 := m2max +. 2.0;
          let lstep = ref 1 in
          let y = ref 0.0 in
          let bbrk = ref false in
          while not !bbrk do
            incr lstep;
            m2 := !m2 -. 1.0;
            let m3 = -.m1 -. !m2 in
            oldfac := !newfac;
            let a1s =
              (l2 -. !m2 +. 2.0) *. (l2 +. !m2 -. 1.0) *. (l3 +. m3 +. 2.0)
              *. (l3 -. m3 -. 1.0)
            in
            newfac := sqrt a1s;
            let dv =
              ((l1 +. l2 +. l3 +. 1.0) *. (l2 +. l3 -. l1))
              -. ((l2 -. !m2 +. 1.0) *. (l3 +. m3 +. 1.0))
              -. ((l2 +. !m2 -. 1.0) *. (l3 -. m3 -. 1.0))
            in
            c1 := -.dv /. !newfac;
            if !lstep > 2 then begin
              let c2 = -. !oldfac /. !newfac in
              y :=
                (!c1 *. thrcof.(nfinp2 - !lstep - 1))
                +. (c2 *. thrcof.(nfinp3 - !lstep - 1));
              if !lstep = nstep2 then bbrk := true
              else begin
                thrcof.(nfin + 1 - !lstep - 1) <- !y;
                sumbac := !sum2;
                sum2 := !sum2 +. (!y *. !y);
                if Float.abs !y >= srhuge then begin
                  for i = 1 to !lstep do
                    let ix = nfin - i + 1 in
                    if Float.abs thrcof.(ix - 1) < srtiny then thrcof.(ix - 1) <- 0.0;
                    thrcof.(ix - 1) <- thrcof.(ix - 1) /. srhuge
                  done;
                  sum2 := !sum2 /. hugee;
                  sumbac := !sumbac /. hugee
                end
              end
            end
            else begin
              y := srtiny *. !c1;
              thrcof.(nfin - 2) <- !y;
              if !lstep = nstep2 then bbrk := true
              else begin
                sumbac := !sum2;
                sum2 := !sum2 +. (!y *. !y)
              end
            end
          done;
          let y3 = !y in
          let y2 = thrcof.(nfinp2 - !lstep - 1) in
          let y1 = thrcof.(nfinp3 - !lstep - 1) in
          let ratio =
            ((x1 *. y1) +. (x2 *. y2) +. (x3 *. y3))
            /. ((x1 *. x1) +. (x2 *. x2) +. (x3 *. x3))
          in
          let nlim = nfin - nstep2 + 1 in
          if Float.abs ratio < 1.0 then begin
            let nlim = nlim + 1 in
            let ratio = 1.0 /. ratio in
            for nn = nlim to nfin do
              thrcof.(nn - 1) <- ratio *. thrcof.(nn - 1)
            done;
            sumuni := !sumfor +. (ratio *. ratio *. !sumbac)
          end
          else begin
            for nn = 1 to nlim do
              thrcof.(nn - 1) <- ratio *. thrcof.(nn - 1)
            done;
            sumuni := (ratio *. ratio *. !sumfor) +. !sumbac
          end;
          brk := true
        end
      end
    done;
    let cnorm = ref (1.0 /. sqrt ((l1 +. l1 +. 1.0) *. !sumuni)) in
    let sign1 = Float.copy_sign 1.0 thrcof.(nfin - 1) in
    let sign2 = parity (int_of_float (Float.abs (l2 -. l3 -. m1) +. eps)) in
    if sign1 *. sign2 <= 0.0 then cnorm := -. !cnorm;
    if Float.abs !cnorm < 1.0 then begin
      let thresh = tinyy /. Float.abs !cnorm in
      for nn = 1 to nfin do
        if Float.abs thrcof.(nn - 1) < thresh then thrcof.(nn - 1) <- 0.0;
        thrcof.(nn - 1) <- !cnorm *. thrcof.(nn - 1)
      done
    end
    else
      for nn = 1 to nfin do
        thrcof.(nn - 1) <- !cnorm *. thrcof.(nn - 1)
      done;
    (m2min, m2max)
  end

let drc6j l2 l3 l4 l5 l6 sixcof ndim =
  if ndim < 1 then invalid_arg "drc6j: ndim < 1";
  if Array.length sixcof < ndim then
    invalid_arg "drc6j: sixcof is shorter than ndim";
  let eps = wigner_eps in
  let hugee = sqrt (huge_dp /. 20.0) in
  let srhuge = sqrt hugee in
  let tinyy = 1.0 /. hugee in
  let srtiny = 1.0 /. srhuge in
  if Float.rem (l2 +. l3 +. l5 +. l6 +. eps) 1.0 >= eps +. eps
     || Float.rem (l4 +. l2 +. l6 +. eps) 1.0 >= eps +. eps
  then invalid_arg "drc6j: l2+l3+l5+l6 or l4+l2+l6 not integer";
  if l4 +. l2 -. l6 < 0.0 || l4 -. l2 +. l6 < 0.0 || -.l4 +. l2 +. l6 < 0.0 then
    invalid_arg "drc6j: l4, l2, l6 triangular condition not satisfied";
  if l4 -. l5 +. l3 < 0.0 || l4 +. l5 -. l3 < 0.0 || -.l4 +. l5 +. l3 < 0.0 then
    invalid_arg "drc6j: l4, l5, l3 triangular condition not satisfied";
  let l1min = Float.max (Float.abs (l2 -. l3)) (Float.abs (l5 -. l6)) in
  let l1max = Float.min (l2 +. l3) (l5 +. l6) in
  if Float.rem (l1max -. l1min +. eps) 1.0 >= eps +. eps then
    invalid_arg "drc6j: l1max-l1min not integer";
  if not (l1min < l1max -. eps) then begin
    if not (l1min < l1max +. eps) then invalid_arg "drc6j: l1min > l1max";
    sixcof.(0) <-
      parity (int_of_float (l2 +. l3 +. l5 +. l6 +. eps))
      /. sqrt ((l1min +. l1min +. 1.0) *. (l4 +. l4 +. 1.0));
    (l1min, l1max)
  end
  else begin
    let nfin = int_of_float (l1max -. l1min +. 1.0 +. eps) in
    if ndim < nfin then
      invalid_arg "drc6j: result array for the 6j coefficients too small";
    let l1 = ref l1min in
    let newfac = ref 0.0 and oldfac = ref 0.0 in
    let c1 = ref 0.0 and c1old = ref 0.0 and denom = ref 0.0 in
    sixcof.(0) <- srtiny;
    let sum1 = ref ((!l1 +. !l1 +. 1.0) *. tinyy) in
    let sum2 = ref 0.0 and sumfor = ref 0.0 and sumbac = ref 0.0 in
    let sumuni = ref 0.0 in
    let lstep = ref 1 in
    let brk = ref false in
    while not !brk do
      incr lstep;
      l1 := !l1 +. 1.0;
      oldfac := !newfac;
      let a1 =
        (!l1 +. l2 +. l3 +. 1.0) *. (!l1 -. l2 +. l3) *. (!l1 +. l2 -. l3)
        *. (-. !l1 +. l2 +. l3 +. 1.0)
      in
      let a2 =
        (!l1 +. l5 +. l6 +. 1.0) *. (!l1 -. l5 +. l6) *. (!l1 +. l5 -. l6)
        *. (-. !l1 +. l5 +. l6 +. 1.0)
      in
      newfac := sqrt (a1 *. a2);
      if !l1 < 1.0 +. eps then
        c1 :=
          -2.0 *. ((l2 *. (l2 +. 1.0)) +. (l5 *. (l5 +. 1.0)) -. (l4 *. (l4 +. 1.0)))
          /. !newfac
      else begin
        let dv =
          (2.0
          *. ((l2 *. (l2 +. 1.0) *. l5 *. (l5 +. 1.0))
             +. (l3 *. (l3 +. 1.0) *. l6 *. (l6 +. 1.0))
             -. (!l1 *. (!l1 -. 1.0) *. l4 *. (l4 +. 1.0))))
          -. (((l2 *. (l2 +. 1.0)) +. (l3 *. (l3 +. 1.0)) -. (!l1 *. (!l1 -. 1.0)))
             *. ((l5 *. (l5 +. 1.0)) +. (l6 *. (l6 +. 1.0))
                -. (!l1 *. (!l1 -. 1.0))))
        in
        denom := (!l1 -. 1.0) *. !newfac;
        if !lstep > 2 then c1old := Float.abs !c1;
        c1 := -.(!l1 +. !l1 -. 1.0) *. dv /. !denom
      end;
      if !lstep <= 2 then begin
        let x = srtiny *. !c1 in
        sixcof.(1) <- x;
        sum1 := !sum1 +. (tinyy *. (!l1 +. !l1 +. 1.0) *. !c1 *. !c1);
        if !lstep = nfin then begin
          sumuni := !sum1;
          brk := true
        end
      end
      else begin
        let c2 = -. !l1 *. !oldfac /. !denom in
        let x = ref ((!c1 *. sixcof.(!lstep - 2)) +. (c2 *. sixcof.(!lstep - 3))) in
        sixcof.(!lstep - 1) <- !x;
        sumfor := !sum1;
        sum1 := !sum1 +. ((!l1 +. !l1 +. 1.0) *. !x *. !x);
        let cycle = ref false in
        if !lstep <> nfin then begin
          if Float.abs !x >= srhuge then begin
            for i = 1 to !lstep do
              if Float.abs sixcof.(i - 1) < srtiny then sixcof.(i - 1) <- 0.0;
              sixcof.(i - 1) <- sixcof.(i - 1) /. srhuge
            done;
            sum1 := !sum1 /. hugee;
            sumfor := !sumfor /. hugee;
            x := !x /. srhuge
          end;
          if !c1old > Float.abs !c1 then cycle := true
        end;
        if not !cycle then begin
          let x1 = !x in
          let x2 = sixcof.(!lstep - 2) in
          let x3 = sixcof.(!lstep - 3) in
          let nfinp2 = nfin + 2 in
          let nfinp3 = nfin + 3 in
          let nstep2 = nfin - !lstep + 3 in
          l1 := l1max;
          sixcof.(nfin - 1) <- srtiny;
          sum2 := (!l1 +. !l1 +. 1.0) *. tinyy;
          l1 := !l1 +. 2.0;
          let lstep = ref 1 in
          let y = ref 0.0 in
          let bbrk = ref false in
          while not !bbrk do
            incr lstep;
            l1 := !l1 -. 1.0;
            oldfac := !newfac;
            let a1s =
              (!l1 +. l2 +. l3) *. (!l1 -. l2 +. l3 -. 1.0)
              *. (!l1 +. l2 -. l3 -. 1.0)
              *. (-. !l1 +. l2 +. l3 +. 2.0)
            in
            let a2s =
              (!l1 +. l5 +. l6) *. (!l1 -. l5 +. l6 -. 1.0)
              *. (!l1 +. l5 -. l6 -. 1.0)
              *. (-. !l1 +. l5 +. l6 +. 2.0)
            in
            newfac := sqrt (a1s *. a2s);
            let dv =
              (2.0
              *. ((l2 *. (l2 +. 1.0) *. l5 *. (l5 +. 1.0))
                 +. (l3 *. (l3 +. 1.0) *. l6 *. (l6 +. 1.0))
                 -. (!l1 *. (!l1 -. 1.0) *. l4 *. (l4 +. 1.0))))
              -. (((l2 *. (l2 +. 1.0)) +. (l3 *. (l3 +. 1.0))
                  -. (!l1 *. (!l1 -. 1.0)))
                 *. ((l5 *. (l5 +. 1.0)) +. (l6 *. (l6 +. 1.0))
                    -. (!l1 *. (!l1 -. 1.0))))
            in
            denom := !l1 *. !newfac;
            c1 := -.(!l1 +. !l1 -. 1.0) *. dv /. !denom;
            if !lstep > 2 then begin
              let c2 = -.(!l1 -. 1.0) *. !oldfac /. !denom in
              y :=
                (!c1 *. sixcof.(nfinp2 - !lstep - 1))
                +. (c2 *. sixcof.(nfinp3 - !lstep - 1));
              if !lstep = nstep2 then bbrk := true
              else begin
                sixcof.(nfin + 1 - !lstep - 1) <- !y;
                sumbac := !sum2;
                sum2 := !sum2 +. ((!l1 +. !l1 -. 3.0) *. !y *. !y);
                if Float.abs !y >= srhuge then begin
                  for i = 1 to !lstep do
                    let ix = nfin - i + 1 in
                    if Float.abs sixcof.(ix - 1) < srtiny then sixcof.(ix - 1) <- 0.0;
                    sixcof.(ix - 1) <- sixcof.(ix - 1) /. srhuge
                  done;
                  sumbac := !sumbac /. hugee;
                  sum2 := !sum2 /. hugee
                end
              end
            end
            else begin
              y := srtiny *. !c1;
              sixcof.(nfin - 2) <- !y;
              if !lstep = nstep2 then bbrk := true
              else begin
                sumbac := !sum2;
                sum2 := !sum2 +. ((!l1 +. !l1 -. 3.0) *. !c1 *. !c1 *. tinyy)
              end
            end
          done;
          let y3 = !y in
          let y2 = sixcof.(nfinp2 - !lstep - 1) in
          let y1 = sixcof.(nfinp3 - !lstep - 1) in
          let ratio =
            ((x1 *. y1) +. (x2 *. y2) +. (x3 *. y3))
            /. ((x1 *. x1) +. (x2 *. x2) +. (x3 *. x3))
          in
          let nlim = nfin - nstep2 + 1 in
          if Float.abs ratio < 1.0 then begin
            let nlim = nlim + 1 in
            let ratio = 1.0 /. ratio in
            for nn = nlim to nfin do
              sixcof.(nn - 1) <- ratio *. sixcof.(nn - 1)
            done;
            sumuni := !sumfor +. (ratio *. ratio *. !sumbac)
          end
          else begin
            for nn = 1 to nlim do
              sixcof.(nn - 1) <- ratio *. sixcof.(nn - 1)
            done;
            sumuni := (ratio *. ratio *. !sumfor) +. !sumbac
          end;
          brk := true
        end
      end
    done;
    let cnorm = ref (1.0 /. sqrt ((l4 +. l4 +. 1.0) *. !sumuni)) in
    let sign1 = Float.copy_sign 1.0 sixcof.(nfin - 1) in
    let sign2 = parity (int_of_float (l2 +. l3 +. l5 +. l6 +. eps)) in
    if sign1 *. sign2 <= 0.0 then cnorm := -. !cnorm;
    if Float.abs !cnorm < 1.0 then begin
      let thresh = tinyy /. Float.abs !cnorm in
      for nn = 1 to nfin do
        if Float.abs sixcof.(nn - 1) < thresh then sixcof.(nn - 1) <- 0.0;
        sixcof.(nn - 1) <- !cnorm *. sixcof.(nn - 1)
      done
    end
    else
      for nn = 1 to nfin do
        sixcof.(nn - 1) <- !cnorm *. sixcof.(nn - 1)
      done;
    (l1min, l1max)
  end
