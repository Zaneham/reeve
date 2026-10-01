(* OCaml side of the gamma differential. It runs the sweep that fdrv.f90 and
   cases.f90 hold, through lib/gamma.ml, and writes the same lines so that the
   two outputs can be compared with cmp. A routine that raises where the
   Fortran returned prints RAISE, which cannot match any value line. *)

open Reeve

let em tag idx v =
  Printf.printf "%-12s %7d %016LX\n" tag idx (Int64.bits_of_float v)

let emi tag idx k = Printf.printf "%-12s %7d %16d\n" tag idx k
let emr tag idx = Printf.printf "%-12s %7d %16s\n" tag idx "RAISE"
let em1 tag idx f = try em tag idx (f ()) with Invalid_argument _ -> emr tag idx

(* ---- The sweep, literal for literal with cases.f90 ---- *)

let lgmcx =
  [| 1.0e3; 2.0e3; 1.0e4; 1.0e5; 1.0e6; 1.0e7; 9.0e7; 9.4e7; 9.48e7; 9.49e7;
     9.5e7; 1.0e8; 1.0e9; 1.0e20; 1.0e100; 1.0e200; 1.0e300; 1.0e305; 1.2e306;
     1.3e306; 1.0e307; 1.0e308 |]

let binn =
  [| 41; 41; 45; 50; 60; 100; 170; 200; 300; 500; 1000; 1000; 1000; 1000; 1000;
     1000; 100000; 100000; 100000; 100000 |]

let binm =
  [| 20; 21; 22; 25; 30; 50; 85; 100; 150; 250; 0; 1; 25; 100; 250; 500; 2; 20;
     21; 40 |]

let pcha =
  [| -3.0; -5.0; -5.0; -5.0; -10.0; -10.0; -15.0; -15.0; -15.0; -15.0; -15.0;
     -20.0; -21.0; -30.0; -30.0; -50.0; -50.0; -100.0; -170.0; -5.0; 0.0; 0.0;
     0.0; -1.0; -1.0; -2.0; -170.0 |]

let pchx =
  [| 0.0; -1.0; -3.0; -5.0; -4.0; 3.0; -3.0; -10.0; 3.0; 10.0; 20.0; -10.0;
     -1.0; -5.0; 5.0; -20.0; 10.0; -30.0; -1.0; 10.0; 0.0; -3.0; 5.0; -1.0;
     0.0; -1.0; 0.0 |]

let poa =
  [| -170.5; -100.5; -60.25; -40.5; -30.5; -25.75; -21.5; -20.5; -16.5; -15.0;
     -13.5; -10.5; -5.5; -2.5; -0.5; -0.25; -0.001; 0.0; 0.001; 0.25; 0.5; 1.0;
     1.5; 2.0; 3.0; 5.5; 10.5; 15.0; 20.5; 21.5; 30.5; 50.5; 100.5; 170.5 |]

let pox =
  [| -30.5; -21.0; -20.5; -20.0; -10.0; -7.5; -5.0; -3.0; -1.0; -0.5; -0.25;
     -0.001; 0.0; 0.001; 0.25; 0.5; 1.0; 2.0; 3.0; 5.0; 7.5; 10.0; 20.0; 20.5;
     21.0; 30.5 |]

let p1a =
  [| -100.5; -50.5; -30.5; -20.5; -10.5; -5.5; -2.5; -1.5; -0.75; -0.5; -0.25;
     0.25; 0.5; 1.0; 2.0; 5.0; 10.0; 11.0; 15.0; 20.0; 50.0; 100.0; 1000.0 |]

let p1x =
  [| 0.0; 1.0e-8; -1.0e-8; 1.0e-5; -1.0e-5; 0.001; -0.001; 0.01; -0.01; 0.02;
     -0.02; 0.05; -0.05; 0.1; -0.1; 0.25; -0.25; 0.5; -0.5; 1.0; -1.0; 2.0;
     -2.0; 5.0; -5.0 |]

let btp =
  [| 0.001; 0.01; 0.1; 0.25; 0.5; 1.0; 1.5; 2.0; 3.0; 5.0; 9.0; 9.5; 9.9; 10.0;
     10.5; 15.0; 20.0; 50.0; 100.0; 500.0; 1000.0; 1.0e4; 1.0e5 |]

let bix =
  [| 0.0; 1.0e-300; 1.0e-20; 1.0e-8; 1.0e-4; 0.001; 0.01; 0.1; 0.19; 0.2; 0.3;
     0.5; 0.7; 0.79; 0.8; 0.9; 0.99; 0.999; 1.0 |]

let biq =
  [| 0.001; 0.1; 0.5; 1.0; 1.5; 2.0; 3.0; 5.0; 9.5; 10.0; 20.0; 50.0; 200.0;
     1000.0 |]

let gia =
  [| 0.001; 0.01; 0.1; 0.5; 1.0; 1.5; 2.0; 3.0; 5.0; 10.0; 20.0; 50.0; 100.0 |]

let gix =
  [| 0.0; 1.0e-10; 1.0e-5; 0.001; 0.01; 0.1; 0.5; 0.9; 1.0; 1.5; 2.0; 5.0;
     10.0; 20.0; 50.0; 100.0; 200.0 |]

let gca =
  [| -20.5; -10.0; -5.001; -5.0; -4.999; -2.5; -1.001; -1.0; -0.999; -0.5;
     -0.001; 0.0; 0.001; 0.5; 1.0; 2.5; 5.0; 10.0; 20.5; 50.0 |]

let gcx =
  [| 1.0e-8; 1.0e-4; 0.001; 0.01; 0.1; 0.5; 0.9; 1.0; 1.5; 2.0; 5.0; 10.0;
     20.0; 50.0; 100.0 |]

let gta =
  [| -20.5; -10.0; -5.0; -4.999; -2.5; -1.0; -0.999; -0.5; -0.001; 0.0; 0.001;
     0.5; 1.0; 2.5; 5.0; 10.0; 20.5; 50.0 |]

let gtx =
  [| 0.0; 1.0e-8; 0.001; 0.01; 0.1; 0.5; 0.9; 1.0; 1.5; 2.0; 5.0; 10.0; 20.0;
     50.0; 100.0 |]

let mica =
  [| 0.0; -0.1; -0.4; -0.5; -0.6; -1.0; -1.5; -2.0; -2.3; -3.0; -5.0; -10.0;
     -20.5 |]

let micx = [| 1.0e-10; 1.0e-5; 0.001; 0.01; 0.1; 0.5; 0.9; 1.0 |]

let mita =
  [| -10.5; -5.5; -2.5; -1.5; -0.999; -0.5; -0.25; 0.0; 0.25; 0.5; 1.0; 2.5;
     5.0; 10.0; 20.5 |]

let mitg =
  [| -12.795895333554364; -2.8130840817693166; 0.86004701537648121;
     1.2655121234846449; 6.9071788853838525; 0.57236494292470042;
     0.20328095143129499; 0.0; -0.098271836421812697; -0.12078223763524543;
     0.0; 1.2009736023470738; 4.7874917427820467; 15.104412573075514;
     43.851925860675159 |]

let mits =
  [| 1.0; -1.0; 1.0; -1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0;
     1.0 |]

let mitx = [| 1.0e-8; 1.0e-4; 0.001; 0.01; 0.1; 0.25; 0.5; 0.9; 1.0 |]

let gica =
  [| -20.5; -10.0; -5.0; -1.0; 0.0; 0.5; 1.0; 2.5; 5.0; 10.0; 20.0; 50.0 |]

let gicx = [| 1.0; 1.5; 2.0; 5.0; 10.0; 20.0; 50.0; 100.0; 500.0 |]
let gita = [| 0.5; 1.0; 2.0; 5.0; 10.0; 20.0; 50.0; 100.0 |]

let gitg =
  [| -0.12078223763524543; 0.0; 0.69314718055994495; 4.7874917427820467;
     15.104412573075514; 42.335616460753485; 148.47776695177305;
     363.73937555556347 |]

let gitf = [| 0.001; 0.01; 0.05; 0.1; 0.25; 0.5; 0.75; 0.9; 0.99; 1.0 |]

let pfx =
  [| 1.0e-20; 1.0e-18; 1.0e-17; 1.0e-16; 1.0e-10; 1.0e-5; 0.001; 0.01; 0.1;
     0.5; 1.0; 2.0; 5.0; 7.0; 10.0; 15.0; 20.0; 50.0; 100.0; 1000.0; 1.0e5;
     1.0e10 |]

let pfn = [| 0; 1; 2; 5; 10; 20; 50; 100; 200 |]
let pfm = [| 1; 2; 3; 5 |]
let chua = [| 0.5; 1.0; 1.5; 2.0; 3.0; -0.5; -1.5; 0.25 |]
let chub = [| 0.3; 0.5; 0.7; 1.0; 1.2; 1.5; 2.0; 2.5; 3.0; 4.4; -0.5 |]
let chuz = [| 1.0; 2.0; 5.0; 10.0; 20.0; 50.0; 100.0; 1000.0; 1.0e4 |]

let chux =
  [| 0.01; 0.1; 0.5; 0.9; 1.0; 1.1; 2.0; 3.0; 5.0; 10.0; 50.0; 100.0 |]

let jjl2 =
  [| 0.0; 1.0; 1.0; 1.0; 2.0; 2.0; 3.0; 3.0; 5.0; 10.0; 10.0; 10.0; 20.0; 50.0;
     100.0; 200.0; 400.0; 0.5; 0.5; 1.5; 1.5; 10.5; 10.5; 20.5; 50.5; 100.5;
     200.5 |]

let jjl3 =
  [| 0.0; 0.0; 1.0; 1.0; 2.0; 3.0; 3.0; 4.0; 5.0; 10.0; 10.0; 10.0; 30.0; 50.0;
     100.0; 150.0; 400.0; 0.5; 1.5; 2.5; 1.5; 20.5; 10.5; 20.5; 50.5; 100.5;
     200.5 |]

let jjm2 =
  [| 0.0; 0.0; 0.0; 1.0; 2.0; 1.0; 0.0; 1.0; 2.0; 0.0; 5.0; 3.0; 5.0; 10.0;
     0.0; 20.0; 1.0; 0.5; 0.5; 0.5; 1.5; 0.5; 10.5; 0.5; 10.5; 0.5; 0.5 |]

let jjm3 =
  [| 0.0; 0.0; 0.0; -1.0; -2.0; -1.0; 0.0; 1.0; -2.0; 0.0; -5.0; 2.0; -5.0;
     -10.0; 0.0; -20.0; -1.0; -0.5; -0.5; -0.5; -1.5; -0.5; -10.5; 0.5; -10.5;
     -0.5; -0.5 |]

let jml1 =
  [| 1.0; 1.0; 1.0; 2.0; 2.0; 5.0; 10.0; 10.0; 20.0; 50.0; 100.0; 200.0; 300.0;
     0.5; 1.5; 1.5; 10.5; 10.5; 50.5; 100.5 |]

let jml2 =
  [| 1.0; 1.0; 0.0; 3.0; 2.0; 5.0; 10.0; 10.0; 30.0; 50.0; 100.0; 150.0; 300.0;
     0.5; 2.5; 1.5; 10.5; 20.5; 50.5; 100.5 |]

let jml3 =
  [| 1.0; 1.0; 1.0; 4.0; 2.0; 5.0; 10.0; 10.0; 40.0; 50.0; 100.0; 100.0; 300.0;
     1.0; 2.0; 1.0; 20.0; 10.0; 100.0; 200.0 |]

let jmm1 =
  [| 0.0; 1.0; 0.0; 1.0; 0.0; 2.0; 0.0; 5.0; 5.0; 10.0; 0.0; 20.0; 0.0; 0.5;
     0.5; 1.5; 0.5; 10.5; 10.5; 0.5 |]

let sjl2 =
  [| 0.0; 1.0; 1.0; 1.0; 2.0; 3.0; 5.0; 10.0; 20.0; 50.0; 100.0; 200.0; 0.5;
     1.5; 1.5; 10.5; 50.5 |]

let sjl3 =
  [| 0.0; 1.0; 1.0; 1.0; 2.0; 4.0; 5.0; 10.0; 30.0; 50.0; 100.0; 200.0; 0.5;
     1.5; 2.5; 10.5; 50.5 |]

let sjl4 =
  [| 0.0; 1.0; 2.0; 0.0; 2.0; 5.0; 5.0; 10.0; 40.0; 50.0; 100.0; 200.0; 1.0;
     1.0; 2.0; 10.0; 50.0 |]

let sjl5 =
  [| 0.0; 1.0; 1.0; 1.0; 2.0; 6.0; 5.0; 10.0; 50.0; 50.0; 100.0; 200.0; 0.5;
     1.5; 2.5; 10.5; 50.5 |]

let sjl6 =
  [| 0.0; 1.0; 1.0; 1.0; 2.0; 7.0; 5.0; 10.0; 60.0; 50.0; 100.0; 200.0; 0.5;
     1.5; 1.5; 10.5; 50.5 |]

(* ---- The three skip rules, mirroring cases.f90 ---- *)

let pochbad a x =
  let ax = a +. x in
  let axint = ax <= 0.0 && Float.trunc ax = ax in
  let aint_ = a <= 0.0 && Float.trunc a = a in
  if axint && not aint_ then true
  else if axint && aint_ && x <> 0.0 then
    Float.min ax a < -20.0 && (a > -9.0 || ax > -9.0)
  else false

let psifn_elim =
  2.302
  *. ((float_of_int (min (-Mach.min_exp_dp) Mach.max_exp_dp)
       *. Mach.log10_radix_dp)
      -. 3.0)

let psifnbad x n m =
  let t = float_of_int (n + m) *. log x in
  Float.abs t > psifn_elim && t <= 0.0

let chuzbad a b =
  let bp = 1.0 +. a -. b in
  bp < 0.0 && Float.trunc bp = bp

let chubad a b x =
  if
    Float.max (Float.abs a) 1.0 *. Float.max (Float.abs (1.0 +. a -. b)) 1.0
    < 0.99 *. Float.abs x
  then chuzbad a b
  else if Float.abs (1.0 +. a -. b) < sqrt Mach.eps_2_dp then true
  else begin
    let aintb =
      if b >= 0.0 then Float.trunc (b +. 0.5) else Float.trunc (b -. 0.5)
    in
    let beps = b -. aintb in
    let n = int_of_float aintb in
    let istrt = if n < 1 then 1 - n else 0 in
    let xi = float_of_int istrt in
    let nonposint v = v <= 0.0 && Float.trunc v = v in
    let bad1 =
      n < 1 && nonposint (1.0 -. b) && not (nonposint (1.0 +. a -. b))
    in
    let bad2 = nonposint (a +. xi -. beps) && not (nonposint a) in
    let bad3 = beps = 0.0 && nonposint (a +. xi) in
    bad1 || bad2 || bad3
  end

(* ---- The sweep ---- *)

let () =
  (* The derived constants, which gfortran folds at compile time and this side
     works out at run time from libm. *)
  em "const" 1 (log Mach.tiny_dp);
  em "const" 2 (log Mach.huge_dp);
  em "const" 3 Mach.log10_radix_dp;
  em "const" 4 Mach.eps_dp;
  em "const" 5 Mach.eps_2_dp;
  em "const" 6 Mach.tiny_dp;
  em "const" 7 Mach.huge_dp;
  em "const" 8 (1.0 /. sqrt Mach.eps_2_dp);
  em "const" 9
    (exp
       (Float.min
          (log (Mach.huge_dp /. 12.0))
          (-.log (12.0 *. Mach.tiny_dp))));
  em "const" 10 (1.0 /. sqrt (24.0 *. Mach.tiny_dp));
  em "const" 11 (log Mach.eps_2_dp);
  em "const" 12 (log Mach.huge_dp -. 0.0001);
  em "const" 13 (0.9 /. Mach.eps_2_dp);
  em "const" 14 (-.log Mach.eps_2_dp);
  em "const" 15 (0.5 *. Mach.eps_2_dp);
  em "const" 16 (sqrt (Mach.huge_dp /. 20.0));
  em "const" 17 (Float.max Mach.eps_dp 1.0e-18);
  em "const" 18 psifn_elim;
  em "const" 19 (Float.max (Mach.eps_dp *. 0.5) 0.5e-18);
  em "const" 20 (Float.min (Mach.log10_radix_dp *. float_of_int Mach.digits_dp) 18.06);

  (* dgamlm *)
  (try
     let xmin, xmax = Gamma.dgamlm () in
     em "dgamlm:min" 1 xmin;
     em "dgamlm:max" 1 xmax
   with Invalid_argument _ ->
     emr "dgamlm:min" 1;
     emr "dgamlm:max" 1);

  (* dgamr and dlgams over x = i/4 from -30 to 170 *)
  let nc = ref 0 in
  for i = -120 to 680 do
    let x = float_of_int i *. 0.25 in
    if not (x <= 0.0 && Float.trunc x = x) then begin
      incr nc;
      em "dgamr:x" !nc x;
      em1 "dgamr" !nc (fun () -> Gamma.dgamr x);
      match Gamma.dlgams x with
      | dl, sg ->
        em "dlgams:l" !nc dl;
        em "dlgams:s" !nc sg
      | exception Invalid_argument _ ->
        emr "dlgams:l" !nc;
        emr "dlgams:s" !nc
    end
  done;

  (* d9lgmc *)
  nc := 0;
  for i = 0 to 400 do
    let x = 10.0 +. (float_of_int i *. 0.5) in
    incr nc;
    em "d9lgmc:x" !nc x;
    em1 "d9lgmc" !nc (fun () -> Gamma.d9lgmc x)
  done;
  Array.iter
    (fun x ->
      incr nc;
      em "d9lgmc:x" !nc x;
      em1 "d9lgmc" !nc (fun () -> Gamma.d9lgmc x))
    lgmcx;

  (* dfac *)
  for n = 0 to 170 do
    emi "dfac:n" (n + 1) n;
    em1 "dfac" (n + 1) (fun () -> Gamma.dfac n)
  done;

  (* dbinom *)
  nc := 0;
  for n = 0 to 40 do
    for m = 0 to n do
      incr nc;
      emi "dbinom:n" !nc n;
      emi "dbinom:m" !nc m;
      em1 "dbinom" !nc (fun () -> Gamma.dbinom n m)
    done
  done;
  Array.iteri
    (fun i n ->
      let m = binm.(i) in
      incr nc;
      emi "dbinom:n" !nc n;
      emi "dbinom:m" !nc m;
      em1 "dbinom" !nc (fun () -> Gamma.dbinom n m))
    binn;

  (* dpoch, the a+x non-positive integer branch first *)
  nc := 0;
  Array.iteri
    (fun i a ->
      let x = pchx.(i) in
      incr nc;
      em "dpoch:a" !nc a;
      em "dpoch:x" !nc x;
      em1 "dpoch" !nc (fun () -> Gamma.dpoch a x))
    pcha;
  Array.iter
    (fun a ->
      Array.iter
        (fun x ->
          if not (pochbad a x) then begin
            incr nc;
            em "dpoch:a" !nc a;
            em "dpoch:x" !nc x;
            em1 "dpoch" !nc (fun () -> Gamma.dpoch a x)
          end)
        pox)
    poa;

  (* dpoch1 *)
  nc := 0;
  Array.iter
    (fun a ->
      Array.iter
        (fun x ->
          if not (pochbad a x) then begin
            incr nc;
            em "dpoch1:a" !nc a;
            em "dpoch1:x" !nc x;
            em1 "dpoch1" !nc (fun () -> Gamma.dpoch1 a x)
          end)
        p1x)
    p1a;

  (* dlbeta and dbeta *)
  nc := 0;
  Array.iter
    (fun p ->
      Array.iter
        (fun q ->
          incr nc;
          em "dlbeta:a" !nc p;
          em "dlbeta:b" !nc q;
          em1 "dlbeta" !nc (fun () -> Gamma.dlbeta p q);
          em1 "dbeta" !nc (fun () -> Gamma.dbeta p q))
        btp)
    btp;

  (* dbetai *)
  nc := 0;
  Array.iter
    (fun x ->
      Array.iter
        (fun p ->
          Array.iter
            (fun q ->
              incr nc;
              em "dbetai:x" !nc x;
              em "dbetai:p" !nc p;
              em "dbetai:q" !nc q;
              em1 "dbetai" !nc (fun () -> Gamma.dbetai x p q))
            biq)
        biq)
    bix;

  (* dgami *)
  nc := 0;
  Array.iter
    (fun a ->
      Array.iter
        (fun x ->
          incr nc;
          em "dgami:a" !nc a;
          em "dgami:x" !nc x;
          em1 "dgami" !nc (fun () -> Gamma.dgami a x))
        gix)
    gia;

  (* dgamic *)
  nc := 0;
  Array.iter
    (fun a ->
      Array.iter
        (fun x ->
          incr nc;
          em "dgamic:a" !nc a;
          em "dgamic:x" !nc x;
          em1 "dgamic" !nc (fun () -> Gamma.dgamic a x))
        gcx)
    gca;
  Array.iter
    (fun a ->
      if a > 0.0 then begin
        incr nc;
        em "dgamic:a" !nc a;
        em "dgamic:x" !nc 0.0;
        em1 "dgamic" !nc (fun () -> Gamma.dgamic a 0.0)
      end)
    gca;

  (* dgamit *)
  nc := 0;
  Array.iter
    (fun a ->
      Array.iter
        (fun x ->
          incr nc;
          em "dgamit:a" !nc a;
          em "dgamit:x" !nc x;
          em1 "dgamit" !nc (fun () -> Gamma.dgamit a x))
        gtx)
    gta;

  (* d9gmic, called directly *)
  nc := 0;
  Array.iter
    (fun a ->
      Array.iter
        (fun x ->
          incr nc;
          em "d9gmic:a" !nc a;
          em "d9gmic:x" !nc x;
          em "d9gmic:l" !nc (log x);
          em1 "d9gmic" !nc (fun () -> Gamma.d9gmic a x (log x)))
        micx)
    mica;

  (* d9gmit, called directly with log|Gamma(a+1)| and its sign handed in as
     literals, so both sides get the same input *)
  nc := 0;
  Array.iteri
    (fun i a ->
      Array.iter
        (fun x ->
          incr nc;
          em "d9gmit:a" !nc a;
          em "d9gmit:x" !nc x;
          em "d9gmit:g" !nc mitg.(i);
          em "d9gmit:s" !nc mits.(i);
          em1 "d9gmit" !nc (fun () -> Gamma.d9gmit a x mitg.(i) mits.(i)))
        mitx)
    mita;

  (* d9lgic, called directly over the a < x region its callers use *)
  nc := 0;
  Array.iter
    (fun a ->
      Array.iter
        (fun x ->
          if a < x then begin
            incr nc;
            em "d9lgic:a" !nc a;
            em "d9lgic:x" !nc x;
            em "d9lgic:l" !nc (log x);
            em1 "d9lgic" !nc (fun () -> Gamma.d9lgic a x (log x))
          end)
        gicx)
    gica;

  (* d9lgit, called directly, needs 0 < x <= a *)
  nc := 0;
  Array.iteri
    (fun i a ->
      Array.iter
        (fun f ->
          let x = f *. a in
          if x > 0.0 && x <= a then begin
            incr nc;
            em "d9lgit:a" !nc a;
            em "d9lgit:x" !nc x;
            em "d9lgit:g" !nc gitg.(i);
            em1 "d9lgit" !nc (fun () -> Gamma.d9lgit a x gitg.(i))
          end)
        gitf)
    gita;

  (* dpsixn *)
  for n = 1 to 200 do
    emi "dpsixn:n" n n;
    em1 "dpsixn" n (fun () -> Gamma.dpsixn n)
  done;

  (* dpsifn, ported from src/original/src/dpsifn.f *)
  nc := 0;
  let ans = Array.make 8 0.0 in
  Array.iter
    (fun x ->
      Array.iter
        (fun n ->
          for kode = 1 to 2 do
            Array.iter
              (fun m ->
                if not (psifnbad x n m) then begin
                  incr nc;
                  Array.fill ans 0 8 0.0;
                  em "dpsifn:x" !nc x;
                  emi "dpsifn:n" !nc n;
                  emi "dpsifn:k" !nc kode;
                  emi "dpsifn:m" !nc m;
                  match Gamma.dpsifn x n kode m ans with
                  | nz ->
                    emi "dpsifn:nz" !nc nz;
                    emi "dpsifn:ie" !nc 0;
                    for jj = 1 to m do
                      em "dpsifn:a" ((!nc * 10) + jj) ans.(jj - 1)
                    done
                  | exception Invalid_argument _ ->
                    emr "dpsifn:nz" !nc;
                    emr "dpsifn:ie" !nc;
                    for jj = 1 to m do
                      emr "dpsifn:a" ((!nc * 10) + jj)
                    done
                end)
              pfm
          done)
        pfn)
    pfx;

  (* d9chu, called directly *)
  nc := 0;
  Array.iter
    (fun a ->
      Array.iter
        (fun b ->
          Array.iter
            (fun z ->
              if not (chuzbad a b) then begin
              incr nc;
              em "d9chu:a" !nc a;
              em "d9chu:b" !nc b;
              em "d9chu:z" !nc z;
              em1 "d9chu" !nc (fun () -> Gamma.d9chu a b z) end)
            chuz)
        chub)
    chua;

  (* dchu *)
  nc := 0;
  Array.iter
    (fun a ->
      Array.iter
        (fun b ->
          Array.iter
            (fun x ->
              if not (chubad a b x) then begin
                incr nc;
                em "dchu:a" !nc a;
                em "dchu:b" !nc b;
                em "dchu:x" !nc x;
                em1 "dchu" !nc (fun () -> Gamma.dchu a b x)
              end)
            chux)
        chub)
    chua;

  let wig = Array.make 4100 0.0 in

  (* drc3jj *)
  Array.iteri
    (fun i0 l2 ->
      let i = i0 + 1 in
      let l3 = jjl3.(i0) and m2 = jjm2.(i0) and m3 = jjm3.(i0) in
      let nf =
        int_of_float
          (l2 +. l3
          -. Float.max (Float.abs (l2 -. l3)) (Float.abs (-.m2 -. m3))
          +. 1.0 +. 0.01)
      in
      em "drc3jj:l2" i l2;
      em "drc3jj:l3" i l3;
      em "drc3jj:m2" i m2;
      em "drc3jj:m3" i m3;
      emi "drc3jj:nd" i nf;
      Array.fill wig 0 4100 0.0;
      match Gamma.drc3jj l2 l3 m2 m3 wig nf with
      | lmn, lmx ->
        em "drc3jj:mn" i lmn;
        em "drc3jj:mx" i lmx;
        for k = 1 to nf do
          em "drc3jj" ((i * 1000) + k) wig.(k - 1)
        done
      | exception Invalid_argument _ ->
        emr "drc3jj:mn" i;
        emr "drc3jj:mx" i;
        for k = 1 to nf do
          emr "drc3jj" ((i * 1000) + k)
        done)
    jjl2;

  (* drc3jm *)
  Array.iteri
    (fun i0 l1 ->
      let i = i0 + 1 in
      let l2 = jml2.(i0) and l3 = jml3.(i0) and m1 = jmm1.(i0) in
      let nf =
        int_of_float
          (Float.min l2 (l3 -. m1)
          -. Float.max (-.l2) (-.l3 -. m1)
          +. 1.0 +. 0.01)
      in
      em "drc3jm:l1" i l1;
      em "drc3jm:l2" i l2;
      em "drc3jm:l3" i l3;
      em "drc3jm:m1" i m1;
      emi "drc3jm:nd" i nf;
      Array.fill wig 0 4100 0.0;
      match Gamma.drc3jm l1 l2 l3 m1 wig nf with
      | lmn, lmx ->
        em "drc3jm:mn" i lmn;
        em "drc3jm:mx" i lmx;
        for k = 1 to nf do
          em "drc3jm" ((i * 1000) + k) wig.(k - 1)
        done
      | exception Invalid_argument _ ->
        emr "drc3jm:mn" i;
        emr "drc3jm:mx" i;
        for k = 1 to nf do
          emr "drc3jm" ((i * 1000) + k)
        done)
    jml1;

  (* drc6j *)
  Array.iteri
    (fun i0 l2 ->
      let i = i0 + 1 in
      let l3 = sjl3.(i0)
      and l4 = sjl4.(i0)
      and l5 = sjl5.(i0)
      and l6 = sjl6.(i0) in
      let nf =
        int_of_float
          (Float.min (l2 +. l3) (l5 +. l6)
          -. Float.max (Float.abs (l2 -. l3)) (Float.abs (l5 -. l6))
          +. 1.0 +. 0.01)
      in
      em "drc6j:l2" i l2;
      em "drc6j:l3" i l3;
      em "drc6j:l4" i l4;
      em "drc6j:l5" i l5;
      em "drc6j:l6" i l6;
      emi "drc6j:nd" i nf;
      Array.fill wig 0 4100 0.0;
      match Gamma.drc6j l2 l3 l4 l5 l6 wig nf with
      | lmn, lmx ->
        em "drc6j:mn" i lmn;
        em "drc6j:mx" i lmx;
        for k = 1 to nf do
          em "drc6j" ((i * 1000) + k) wig.(k - 1)
        done
      | exception Invalid_argument _ ->
        emr "drc6j:mn" i;
        emr "drc6j:mx" i;
        for k = 1 to nf do
          emr "drc6j" ((i * 1000) + k)
        done)
    sjl2
