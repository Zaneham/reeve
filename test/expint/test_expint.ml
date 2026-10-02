(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Level 1, level 2 and level 3 checks for the exponential integrals, the
   logarithmic integral, Spence's integral, Dawson's function and the
   elementary kernels. *)

open Reeve.Expint
open Reeve.Expint.Raw

let pi = 3.14159265358979323846
let ln2 = 0.69314718055994530942
let tol = 1.0e-10
let loose_tol = 1.0e-6
let hist_tol = 1.0e-15
let bad = ref 0

let fail what =
  Printf.printf "FAIL %s\n" what;
  bad := !bad + 1

let near what a b = if not (Float.abs (a -. b) < tol) then fail what
let loose what a b = if not (Float.abs (a -. b) < loose_tol) then fail what
let near_t what t a b = if not (Float.abs (a -. b) < t) then fail what

let rel_t what t a b =
  let e =
    if Float.abs b > 1.0e-14 then Float.abs (a -. b) /. Float.abs b
    else Float.abs (a -. b)
  in
  if not (e < t) then fail what

let holds what c = if not c then fail what
let bits what a b = if Int64.bits_of_float a <> Int64.bits_of_float b then fail what

(* ---- Level 1, regression ---- *)

let l1_expint () =
  loose "l1 de1 1" (de1 1.0) 0.21938393439552029;
  loose "l1 de1 2" (de1 2.0) 0.04890051070806112;
  loose "l1 dei 1" (dei 1.0) 1.8951178163559368;
  near_t "l1 de1 x plus dei -x" loose_tol (de1 2.0 +. dei (-2.0)) 0.0;
  holds "l1 de1 decreases for x > 0"
    (de1 1.0 > de1 2.0 && de1 2.0 > 0.0);
  holds "l1 de1 underflows far out" (de1 1.0e4 = 0.0)

let l1_others () =
  near "l1 dli 2" (dli 2.0) 1.0451637801174929;
  near "l1 dspenc 0" (dspenc 0.0) 0.0;
  near "l1 ddaws 0" (ddaws 0.0) 0.0;
  near "l1 dexprl 0" (dexprl 0.0) 1.0;
  near "l1 dlnrel 0" (dlnrel 0.0) 0.0;
  near "l1 dcbrt 27" (dcbrt 27.0) 3.0;
  near "l1 dcbrt -8" (dcbrt (-8.0)) (-2.0);
  near "l1 dsindg 0" (dsindg 0.0) 0.0;
  near "l1 dcosdg 0" (dcosdg 0.0) 1.0;
  near "l1 d9atn1 0" (d9atn1 0.0) (-1.0 /. 3.0);
  near "l1 d9ln2r 0" (d9ln2r 0.0) (1.0 /. 3.0);
  near "l1 d9pak 0.5 3" (d9pak 0.5 3) 4.0

(* ---- Level 2, mathematical identities ---- *)

let l2_continuation () =
  List.iter
    (fun x ->
      bits (Printf.sprintf "l2 de1 is minus dei reflected at %g" x) (de1 x)
        (-.dei (-.x)))
    [ 0.25; 1.0; 2.5; 7.0; 40.0; -0.5; -3.0; -20.0 ];
  List.iter
    (fun x -> bits (Printf.sprintf "l2 dli is dei of log at %g" x) (dli x)
        (dei (log x)))
    [ 0.1; 0.5; 1.5; 2.0; 10.0; 1.0e6 ]

let l2_spence () =
  rel_t "l2 dspenc 1" hist_tol (dspenc 1.0) (pi *. pi /. 6.0);
  rel_t "l2 dspenc 0.5" hist_tol (dspenc 0.5)
    ((pi *. pi /. 12.0) -. (ln2 *. ln2 /. 2.0));
  rel_t "l2 dspenc 2" hist_tol (dspenc 2.0) (pi *. pi /. 4.0);
  rel_t "l2 dspenc -1" hist_tol (dspenc (-1.0)) (-.pi *. pi /. 12.0);
  let x = 0.3 in
  let series =
    let s = ref 0.0 and t = ref 1.0 in
    for k = 1 to 60 do
      t := !t *. x;
      s := !s +. (!t /. float_of_int (k * k))
    done;
    !s
  in
  rel_t "l2 dspenc 0.3 against its series" 1.0e-14 (dspenc x) series

let l2_dawson () =
  List.iter
    (fun x ->
      let h = 1.0e-5 in
      let d = (ddaws (x +. h) -. ddaws (x -. h)) /. (2.0 *. h) in
      near_t
        (Printf.sprintf "l2 ddaws ode at %g" x)
        1.0e-8 d
        (1.0 -. (2.0 *. x *. ddaws x)))
    [ 0.25; 0.5; 1.0; 2.0; 3.5 ];
  List.iter
    (fun x ->
      bits (Printf.sprintf "l2 ddaws odd at %g" x) (ddaws (-.x)) (-.ddaws x))
    [ 0.1; 0.5; 1.0; 2.0; 5.0; 100.0 ];
  near_t "l2 ddaws asymptotic" 1.0e-6 (ddaws 1.0e4) (0.5 /. 1.0e4)

let l2_elementary () =
  List.iter
    (fun x ->
      rel_t
        (Printf.sprintf "l2 dlnrel inverts dexprl at %g" x)
        1.0e-14
        (dlnrel (x *. dexprl x))
        x)
    [ 1.0e-12; 1.0e-9; 1.0e-6; 1.0e-3; 0.01; 0.1; 0.3; -1.0e-9; -1.0e-3; -0.1;
      -0.3 ];
  List.iter
    (fun x ->
      rel_t
        (Printf.sprintf "l2 dexprl is the exp quotient at %g" x)
        1.0e-14 (dexprl x)
        ((exp x -. 1.0) /. x))
    [ 0.6; 1.0; 2.0; 5.0 ];
  List.iter
    (fun x ->
      rel_t
        (Printf.sprintf "l2 d9atn1 reconstructs atan at %g" x)
        1.0e-13
        (x +. (x *. x *. x *. d9atn1 x))
        (atan x))
    [ 1.0e-6; 0.1; 0.5; 0.9; 1.0; 2.0; 100.0; -0.5; -1.5 ];
  List.iter
    (fun x ->
      rel_t
        (Printf.sprintf "l2 d9ln2r reconstructs log1p at %g" x)
        1.0e-14
        (x -. (0.5 *. x *. x) +. (x *. x *. x *. d9ln2r x))
        (dlnrel x))
    [ 1.0e-6; 0.1; 0.5; 0.8; 0.9; 2.0; -0.1; -0.5; -0.6; -0.7 ];
  List.iter
    (fun x ->
      rel_t (Printf.sprintf "l2 dcbrt cubes back at %g" x) 1.0e-14
        (dcbrt x *. dcbrt x *. dcbrt x)
        x)
    [ 1.0e-300; 1.0e-20; 0.5; 2.0; 10.0; 1.0e20; 1.0e300; -1.0e-20; -2.0;
      -1.0e300 ];
  bits "l2 dcbrt of a zero" (dcbrt 0.0) 0.0

let l2_degrees () =
  bits "l2 dsindg 90" (dsindg 90.0) 1.0;
  bits "l2 dsindg 180" (dsindg 180.0) 0.0;
  bits "l2 dsindg 270" (dsindg 270.0) (-1.0);
  bits "l2 dsindg -90" (dsindg (-90.0)) (-1.0);
  bits "l2 dcosdg 0" (dcosdg 0.0) 1.0;
  bits "l2 dcosdg 90" (dcosdg 90.0) 0.0;
  bits "l2 dcosdg 180" (dcosdg 180.0) (-1.0);
  bits "l2 dcosdg 360" (dcosdg 360.0) 1.0;
  List.iter
    (fun x ->
      rel_t
        (Printf.sprintf "l2 degree pythagoras at %g" x)
        1.0e-15
        ((dsindg x *. dsindg x) +. (dcosdg x *. dcosdg x))
        1.0;
      rel_t
        (Printf.sprintf "l2 dsindg against radians at %g" x)
        1.0e-14 (dsindg x)
        (sin (x *. pi /. 180.0)))
    [ 1.0; 17.0; 30.0; 45.0; 60.0; 123.0; 200.0; -37.0 ]

let l2_pack () =
  List.iter
    (fun x ->
      let y, n = d9upak x in
      holds
        (Printf.sprintf "l2 d9upak normalises %g" x)
        (Float.abs y >= 0.5 && Float.abs y < 1.0);
      bits (Printf.sprintf "l2 d9pak inverts d9upak at %g" x) (d9pak y n) x)
    [ 1.0; 0.5; 0.1; 3.0; 1.0e-300; 1.0e300; Float.max_float; -7.25;
      Float.min_float ];
  let y, n = d9upak 0.0 in
  holds "l2 d9upak of a zero" (y = 0.0 && n = 0);
  bits "l2 d9pak scales by a power of two" (d9pak 0.75 10) (0.75 *. 1024.0);
  bits "l2 d9pak underflows to zero rather than denormalising"
    (d9pak 0.5 (-1073)) 0.0

let l2_sequence () =
  let m = 5 in
  List.iter
    (fun (x, n) ->
      let en = Array.make m 0.0 in
      let nz = dexint x n 1 m 1.0e-14 en in
      holds (Printf.sprintf "l2 dexint nz zero at x=%g n=%d" x n) (nz = 0);
      let emx = exp (-.x) in
      for k = 0 to m - 2 do
        let kk = float_of_int (n + k) in
        rel_t
          (Printf.sprintf "l2 dexint recurrence at x=%g n=%d k=%d" x n k)
          1.0e-12
          ((kk *. en.(k + 1)) +. (x *. en.(k)))
          emx
      done)
    [ (0.25, 1); (0.5, 2); (1.0, 1); (1.5, 3); (2.5, 1); (3.0, 2); (5.0, 1);
      (10.0, 4); (20.0, 2); (60.0, 1) ];
  let en = Array.make 1 0.0 in
  let _ = dexint 1.0 1 1 1 1.0e-14 en in
  rel_t "l2 dexint n=1 against de1" 1.0e-13 en.(0) (de1 1.0);
  let en = Array.make 4 0.0 in
  let _ = dexint 0.0 2 1 4 1.0e-14 en in
  Array.iteri
    (fun i v ->
      rel_t
        (Printf.sprintf "l2 dexint at x=0 index %d" i)
        1.0e-15 v
        (1.0 /. float_of_int (i + 1)))
    en;
  let en = Array.make 3 1.0 in
  let nz = dexint 705.0 1 1 3 1.0e-14 en in
  holds "l2 dexint underflow flag" (nz = 3);
  holds "l2 dexint underflow fills zeros"
    (en.(0) = 0.0 && en.(1) = 0.0 && en.(2) = 0.0);
  let en = Array.make 3 0.0 in
  let _ = dexint 4.0 1 2 3 1.0e-14 en in
  let es = Array.make 3 0.0 in
  let _ = dexint 4.0 1 1 3 1.0e-14 es in
  for k = 0 to 2 do
    rel_t
      (Printf.sprintf "l2 dexint kode 2 scales index %d" k)
      1.0e-14 en.(k)
      (es.(k) *. exp 4.0)
  done

(* ---- Level 3, values from the reference Fortran ---- *)

let l3_expint () =
  rel_t "l3 de1 1" hist_tol (de1 1.0) 2.19383934395520258e-01;
  rel_t "l3 de1 2" hist_tol (de1 2.0) 4.89005107080611248e-02;
  rel_t "l3 de1 0.5" hist_tol (de1 0.5) 5.59773594776160843e-01;
  rel_t "l3 de1 10" hist_tol (de1 10.0) 4.15696892968532464e-06;
  rel_t "l3 de1 -1" hist_tol (de1 (-1.0)) (-1.89511781635593679e+00);
  rel_t "l3 de1 -5" hist_tol (de1 (-5.0)) (-4.01852753558031779e+01);
  rel_t "l3 dei 1" hist_tol (dei 1.0) 1.89511781635593679e+00;
  rel_t "l3 dei 2" hist_tol (dei 2.0) 4.95423435600189066e+00;
  rel_t "l3 dei -1" hist_tol (dei (-1.0)) (-2.19383934395520258e-01);
  rel_t "l3 dli 2" hist_tol (dli 2.0) 1.04516378011749289e+00;
  rel_t "l3 dli 10" hist_tol (dli 10.0) 6.16559950478730023e+00;
  rel_t "l3 dli 1e6" hist_tol (dli 1.0e6) 7.86275491594621562e+04;
  rel_t "l3 dli 0.5" hist_tol (dli 0.5) (-3.78671043061087953e-01)

let l3_others () =
  rel_t "l3 dspenc 0.5" hist_tol (dspenc 0.5) 5.82240526465012453e-01;
  rel_t "l3 dspenc 1" hist_tol (dspenc 1.0) 1.64493406684822641e+00;
  rel_t "l3 dspenc 2" hist_tol (dspenc 2.0) 2.46740110027233950e+00;
  rel_t "l3 dspenc -1" hist_tol (dspenc (-1.0)) (-8.22467033424113203e-01);
  rel_t "l3 dspenc 0.25" hist_tol (dspenc 0.25) 2.67652639082732624e-01;
  rel_t "l3 dspenc 4" hist_tol (dspenc 4.0) 2.06130946677731774e+00;
  rel_t "l3 ddaws 0.5" hist_tol (ddaws 0.5) 4.24436383502022285e-01;
  rel_t "l3 ddaws 1" hist_tol (ddaws 1.0) 5.38079506912768402e-01;
  rel_t "l3 ddaws 2" hist_tol (ddaws 2.0) 3.01340388923791946e-01;
  rel_t "l3 ddaws 5" hist_tol (ddaws 5.0) 1.02134074424276841e-01;
  rel_t "l3 ddaws at its maximum" hist_tol
    (ddaws 0.9241388730)
    5.41044224635181759e-01;
  rel_t "l3 dexprl 0.1" hist_tol (dexprl 0.1) 1.05170918075647624e+00;
  rel_t "l3 dexprl -0.25" hist_tol (dexprl (-0.25)) 8.84796867714380486e-01;
  rel_t "l3 dlnrel 0.1" hist_tol (dlnrel 0.1) 9.53101798043248655e-02;
  rel_t "l3 dlnrel -0.25" hist_tol (dlnrel (-0.25)) (-2.87682072451780901e-01);
  rel_t "l3 dcbrt 2" hist_tol (dcbrt 2.0) 1.25992104989487341e+00;
  rel_t "l3 dcbrt -0.001" hist_tol (dcbrt (-0.001)) (-9.99999999999999917e-02);
  rel_t "l3 dsindg 30" hist_tol (dsindg 30.0) 4.99999999999999944e-01;
  rel_t "l3 dcosdg 60" hist_tol (dcosdg 60.0) 5.00000000000000111e-01;
  rel_t "l3 d9atn1 0.5" hist_tol (d9atn1 0.5) (-2.90819127993551085e-01);
  rel_t "l3 d9atn1 -0.75" hist_tol (d9atn1 (-0.75)) (-2.52441816193696267e-01);
  rel_t "l3 d9ln2r 0.5" hist_tol (d9ln2r 0.5) 2.43720864865315051e-01;
  rel_t "l3 d9ln2r -0.5" hist_tol (d9ln2r (-0.5)) 5.45177444479562512e-01

let l3_sequence () =
  let en = Array.make 4 0.0 in
  let _ = dexint 1.0 1 1 4 1.0e-14 en in
  rel_t "l3 dexint 1 n1 k1 1" hist_tol en.(0) 2.19383934395520508e-01;
  rel_t "l3 dexint 1 n1 k1 2" hist_tol en.(1) 1.48495506775921826e-01;
  rel_t "l3 dexint 1 n1 k1 3" hist_tol en.(2) 1.09691967197760254e-01;
  rel_t "l3 dexint 1 n1 k1 4" hist_tol en.(3) 8.60624913245606887e-02;
  let en = Array.make 3 0.0 in
  let _ = dexint 5.0 2 1 3 1.0e-14 en in
  rel_t "l3 dexint 5 n2 k1 1" hist_tol en.(0) 9.96469042708838029e-04;
  rel_t "l3 dexint 5 n2 k1 2" hist_tol en.(1) 8.77800892770638211e-04;
  rel_t "l3 dexint 5 n2 k1 3" hist_tol en.(2) 7.82980845077425208e-04;
  let en = Array.make 3 0.0 in
  let _ = dexint 0.5 3 2 3 1.0e-14 en in
  rel_t "l3 dexint 0.5 n3 k2 1" hist_tol en.(0) 3.65363829060466327e-01;
  rel_t "l3 dexint 0.5 n3 k2 2" hist_tol en.(1) 2.72439361823255621e-01;
  rel_t "l3 dexint 0.5 n3 k2 3" hist_tol en.(2) 2.15945079772093040e-01

(* ---- Domains the Fortran refuses ---- *)

let raises what f =
  match f () with
  | _ -> fail what
  | exception Invalid_argument _ -> ()

let l1_domains () =
  raises "de1 0 raises" (fun () -> de1 0.0);
  raises "dei 0 raises" (fun () -> dei 0.0);
  raises "dli 0 raises" (fun () -> dli 0.0);
  raises "dli 1 raises" (fun () -> dli 1.0);
  raises "dli negative raises" (fun () -> dli (-1.0));
  raises "dlnrel -1 raises" (fun () -> dlnrel (-1.0));
  raises "dlnrel -2 raises" (fun () -> dlnrel (-2.0));
  raises "d9atn1 too big raises" (fun () -> d9atn1 1.0e18);
  raises "d9ln2r too big raises" (fun () -> d9ln2r 1.0e12);
  raises "d9pak overflow raises" (fun () -> d9pak 0.75 1100);
  let en = Array.make 4 0.0 in
  raises "dexint negative x raises" (fun () -> dexint (-1.0) 1 1 4 1.0e-14 en);
  raises "dexint n zero raises" (fun () -> dexint 1.0 0 1 4 1.0e-14 en);
  raises "dexint bad kode raises" (fun () -> dexint 1.0 1 3 4 1.0e-14 en);
  raises "dexint m zero raises" (fun () -> dexint 1.0 1 1 0 1.0e-14 en);
  raises "dexint loose tol raises" (fun () -> dexint 1.0 1 1 4 0.5 en);
  raises "dexint tight tol raises" (fun () -> dexint 1.0 1 1 4 1.0e-20 en);
  raises "dexint x zero with n one raises" (fun () ->
      dexint 0.0 1 1 4 1.0e-14 en)

(* ---- The OCaml surface ---- *)

let surface () =
  let e = sequence ~n:1 ~count:4 1.0 in
  holds "surface sequence length" (Array.length e = 4);
  rel_t "surface sequence against de1" tol e.(0) (de1 1.0);
  rel_t "surface sequence 1 n1 2" hist_tol e.(1) 1.48495506775921826e-01;
  rel_t "surface sequence 1 n1 3" hist_tol e.(2) 1.09691967197760254e-01;
  rel_t "surface sequence 1 n1 4" hist_tol e.(3) 8.60624913245606887e-02;
  let raw = Array.make 3 0.0 in
  let nz = dexint 5.0 2 1 3 1.0e-14 raw in
  holds "surface raw no underflow at 5" (nz = 0);
  let s = sequence ~tol:1.0e-14 ~n:2 ~count:3 5.0 in
  holds "surface agrees with raw bit for bit"
    (Int64.bits_of_float s.(0) = Int64.bits_of_float raw.(0)
    && Int64.bits_of_float s.(1) = Int64.bits_of_float raw.(1)
    && Int64.bits_of_float s.(2) = Int64.bits_of_float raw.(2));
  let u = sequence ~tol:1.0e-14 ~n:1 ~count:3 4.0 in
  let k = sequence ~tol:1.0e-14 ~scaling:Exp_scaled ~n:1 ~count:3 4.0 in
  for i = 0 to 2 do
    rel_t
      (Printf.sprintf "surface exp scaled %d" i)
      1.0e-14 k.(i)
      (u.(i) *. exp 4.0)
  done;
  let z = sequence ~n:1 ~count:3 705.0 in
  holds "surface underflow is all zeros"
    (z.(0) = 0.0 && z.(1) = 0.0 && z.(2) = 0.0);
  let z0 = sequence ~n:2 ~count:4 0.0 in
  for i = 0 to 3 do
    near
      (Printf.sprintf "surface at zero %d" i)
      z0.(i)
      (1.0 /. float_of_int (i + 1))
  done;
  let one = sequence ~n:3 ~count:1 0.25 in
  holds "surface count one" (Array.length one = 1);
  rel_t "surface count one value" tol one.(0)
    (let en = Array.make 1 0.0 in
     ignore (dexint 0.25 3 1 1 1.0e-14 en);
     en.(0));
  raises "surface negative x raises" (fun () -> sequence ~n:1 ~count:4 (-1.0));
  raises "surface n zero raises" (fun () -> sequence ~n:0 ~count:4 1.0);
  raises "surface count zero raises" (fun () -> sequence ~n:1 ~count:0 1.0);
  raises "surface loose tol raises" (fun () ->
      sequence ~tol:0.5 ~n:1 ~count:4 1.0);
  raises "surface tight tol raises" (fun () ->
      sequence ~tol:1.0e-20 ~n:1 ~count:4 1.0);
  raises "surface x zero with n one raises" (fun () ->
      sequence ~n:1 ~count:4 0.0)

let () =
  l1_expint ();
  l1_others ();
  l1_domains ();
  l2_continuation ();
  l2_spence ();
  l2_dawson ();
  l2_elementary ();
  l2_degrees ();
  l2_pack ();
  l2_sequence ();
  l3_expint ();
  l3_others ();
  l3_sequence ();
  surface ();
  if !bad = 0 then print_string "expint: all checks passed\n" else exit 1
