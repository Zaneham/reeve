(* Reeve differential driver for the extended-range Legendre package. It
   writes one input file per DXSET configuration, because the SLATEC DXSET
   refuses to run twice in a process, and writes its own answers to ml.out in
   the same order. *)

open Reeve.Xlegf
open Reeve.Xlegf.Raw

let hex v = Printf.sprintf "%016LX" (Int64.bits_of_float v)
let pi = 4.0 *. atan 1.0
let deg01 = 0.1 *. pi /. 180.0

type cmd =
  | Adj of float * int
  | Rt of float * int
  | Add of float * int * float * int
  | Red of float * int
  | Con of float * int
  | Legf of float * int * int * int * float * int
  | Nrmp of int * int * int * float * int

let line = function
  | Adj (x, ix) -> Printf.sprintf "adj %s %d" (hex x) ix
  | Rt (x, ix) -> Printf.sprintf "rt %s %d" (hex x) ix
  | Add (x, ix, y, iy) -> Printf.sprintf "add %s %d %s %d" (hex x) ix (hex y) iy
  | Red (x, ix) -> Printf.sprintf "red %s %d" (hex x) ix
  | Con (x, ix) -> Printf.sprintf "con %s %d" (hex x) ix
  | Legf (d, nd, m1, m2, th, id) ->
    Printf.sprintf "legf %s %d %d %d %s %d" (hex d) nd m1 m2 (hex th) id
  | Nrmp (nu, m1, m2, da, md) ->
    Printf.sprintf "nrmp %d %d %d %s %d" nu m1 m2 (hex da) md

(* ---- The DXSET configurations ---- *)

let configs =
  [
    (0, 0, 0.0, 0);
    (2, 53, 0.0, 62);
    (2, 53, 0.0, 31);
    (2, 53, 0.0, 20);
    (2, 53, 0.0, 16);
    (2, 53, 0.0, 15);
    (2, 24, 0.0, 24);
    (2, 100, 0.0, 62);
    (2, 53, 1e200, 62);
    (2, 10, 1e40, 62);
    (4, 26, 1e150, 62);
    (8, 17, 1e150, 62);
    (16, 13, 1e150, 62);
    (3, 0, 0.0, 0);
    (2, 53, 0.0, 14);
    (2, 53, 0.0, 64);
    (2, 53, 100.0, 0);
    (2, 600, 0.0, 62);
    (2, 53, 1e-300, 62);
  ]

(* The scaling exponent, the index bound and the radix power, worked out the
   way DXSET works them out. They only choose where the sweep puts its test
   points; nothing here is compared against anything. *)
let derive (irad, nradpl, dzero, nbits) =
  let iradx = if irad = 0 then 2 else irad in
  let nrdplc = if nradpl = 0 then 53 else nradpl in
  let nbitsx = if nbits = 0 then Sys.int_size - 1 else nbits in
  let log2r = match iradx with 2 -> 1 | 4 -> 2 | 8 -> 3 | 16 -> 4 | _ -> 1 in
  let dlg10r = log10 (float_of_int iradx) in
  let l =
    if dzero <> 0.0 then int_of_float (0.5 *. log10 dzero /. dlg10r) - 1
    else min ((1 - -1021) / 2) ((1024 - 1) / 2)
  in
  let l = if l < 4 then 4 else l in
  let nbitsx = if nbitsx < 15 || nbitsx > 62 then 62 else nbitsx in
  let kmax = (1 lsl (nbitsx - 1)) - (2 * l) in
  (log2r, nrdplc, l, kmax)

(* ---- Test points ---- *)

let dedup l =
  let seen = Hashtbl.create 97 in
  List.filter
    (fun v ->
      if Hashtbl.mem seen v then false
      else begin
        Hashtbl.add seen v ();
        true
      end)
    l

let both v = [ v; -.v ]

let principals log2r l =
  let radixl = ldexp 1.0 (log2r * l) in
  let rad2l = radixl *. radixl in
  let ir = 1.0 /. radixl and ir2 = 1.0 /. rad2l in
  0.0 :: -0.0
  :: dedup
       (List.concat_map both
         [
           1.0; Float.pred 1.0; Float.succ 1.0; 0.5; 1.5; 3.25; 123.456;
           radixl; Float.pred radixl; Float.succ radixl; radixl *. 0.5;
           ir; Float.pred ir; Float.succ ir; ir *. 2.0;
           rad2l; Float.pred rad2l; ir2; Float.succ ir2; ir2 *. 2.0;
           1e-300; 1e300; 1e-8; 1e8; min_float; 4.9e-324; max_float;
         ])

let indices l kmax =
  let l2 = 2 * l in
  dedup
    (0
    :: List.concat_map
         (fun k -> [ k; -k ])
         [
           1; l - 1; l; l + 1; l2 - 1; l2; l2 + 1; 3 * l; 6 * l - 1; 6 * l;
           6 * l + 1; 12 * l; kmax / 2; kmax - l2 - 1; kmax - l2; kmax - l2 + 1;
           kmax - 1; kmax; kmax + 1;
         ])

let add_principals log2r l =
  let radixl = ldexp 1.0 (log2r * l) in
  0.0 :: -0.0
  :: dedup
       (List.concat_map both
          [ 1.0; Float.pred 1.0; 0.5; 3.25; radixl *. 0.5; 1.0 /. radixl; 1e-300 ])

let add_indices l kmax =
  let l2 = 2 * l in
  dedup
    (0
    :: List.concat_map
         (fun k -> [ k; -k ])
         [ 1; l; l2; 3 * l; 6 * l; 6 * l + 1; 12 * l; kmax - l2; kmax ])

let con_indices l kmax =
  let l2 = 2 * l in
  dedup
    (0
    :: List.concat_map
         (fun k -> [ k; -k ])
         [
           1; l; l2; 2 * l2; 38 * l2; 100 * l2; 6 * l; 6 * l + 1; kmax - l2;
           kmax;
         ])

(* ---- The Legendre cases ---- *)

let thetas = [ deg01; 0.001; 0.1; 0.5; 1.0; 1.5; pi /. 2.0 ]

let nuwise =
  [
    (0.0, 0, 0); (0.0, 5, 0); (0.0, 0, 1); (0.0, 0, 2); (0.0, 0, 5);
    (1.0, 0, 0); (1.0, 3, 1); (2.0, 4, 2); (10.0, 6, 3); (0.5, 0, 0);
    (0.5, 4, 2); (0.3, 0, 0); (0.3, 3, 1); (10.3, 4, 3); (-0.5, 0, 0);
    (-0.5, 2, 1); (-0.5, 0, 4); (100.0, 0, 10); (100.0, 5, 10);
    (100.0, 0, 100); (2000.4, 5, 2000); (2000.4, 0, 2000); (5000.0, 3, 5000);
    (0.0, 0, 40); (1.0, 0, 30);
  ]

let muwise =
  [
    (0.0, 0, 5); (0.0, 0, 1); (0.0, 2, 7); (1.0, 0, 4); (2.0, 0, 8);
    (0.5, 0, 5); (0.3, 0, 6); (0.3, 2, 8); (-0.5, 0, 3); (10.3, 0, 12);
    (100.0, 0, 20); (100.0, 10, 30); (2000.4, 1995, 2000);
    (5000.0, 4990, 5000); (0.0, 10, 50);
  ]

let legf_bad =
  [
    (1.0, -1, 0, 0, 1.0, 1); (-0.6, 0, 0, 0, 1.0, 1); (1.0, 0, 3, 2, 1.0, 1);
    (1.0, 0, -1, 0, 1.0, 1); (1.0, 0, 0, 0, 0.0, 1);
    (1.0, 0, 0, 0, (pi /. 2.0) +. 0.1, 1); (1.0, 0, 0, 0, -1.0, 1);
    (1.0, 2, 0, 3, 1.0, 1); (0.3, 0, 0, 0, 1.0, 4); (0.3, 2, 0, 0, 1.0, 4);
    (1.0, 0, 0, 0, pi, 2); (0.5, 0, 0, 0, 1.0, 4);
  ]

let nrmp_nus = [ 0; 1; 2; 5; 12; 40; 100; 500; 1000 ]

let nrmp_mus nu =
  dedup
    [ (0, 0); (0, nu / 2); (0, min nu 60); (nu, nu); (1, 3); (3, 4);
      (nu + 1, nu + 2) ]

let nrmp_args =
  [
    (1, 0.0); (1, 0.3); (1, -0.3); (1, 0.9); (1, -0.9); (1, 1.0); (1, -1.0);
    (1, 0.999999); (1, 1e-300); (2, 0.0); (2, deg01); (2, 0.1); (2, 1.0);
    (2, pi /. 2.0); (2, pi); (2, -.pi); (2, 3.0); (2, -0.5);
  ]

let nrmp_bad =
  [
    (-1, 0, 0, 0.3, 1); (1, -1, 0, 0.3, 1); (1, 3, 2, 0.3, 1);
    (5, 0, 0, 1.5, 1); (5, 0, 0, -1.5, 1); (5, 0, 0, 4.0, 2);
    (5, 0, 0, -4.0, 2);
  ]

let cases cfg =
  let log2r, _nrdplc, l, kmax = derive cfg in
  let ps = principals log2r l and ix = indices l kmax in
  let aps = add_principals log2r l and aix = add_indices l kmax in
  let cix = con_indices l kmax in
  let out = ref [] in
  let put c = out := c :: !out in
  List.iter (fun v -> List.iter (fun k -> put (Adj (v, k))) ix) ps;
  List.iter (fun v -> List.iter (fun k -> put (Rt (v, k))) ix) ps;
  List.iter (fun v -> List.iter (fun k -> put (Red (v, k))) ix) ps;
  List.iter (fun v -> List.iter (fun k -> put (Con (v, k))) cix) ps;
  List.iter
    (fun x ->
      List.iter
        (fun i ->
          List.iter
            (fun y -> List.iter (fun j -> put (Add (x, i, y, j))) aix)
            aps)
        aix)
    aps;
  List.iter
    (fun (d, nd, m) ->
      if m + nd + 1 <= 120 then
        List.iter
          (fun th ->
            for id = 1 to 4 do
              put (Legf (d, nd, m, m, th, id))
            done)
          thetas)
    nuwise;
  List.iter
    (fun (d, m1, m2) ->
      if m2 - m1 + 1 <= 120 then
        List.iter
          (fun th ->
            for id = 1 to 4 do
              put (Legf (d, 0, m1, m2, th, id))
            done)
          thetas)
    muwise;
  List.iter
    (fun (d, nd, m1, m2, th, id) -> put (Legf (d, nd, m1, m2, th, id)))
    legf_bad;
  List.iter
    (fun nu ->
      List.iter
        (fun (m1, m2) ->
          if m2 >= m1 && m2 - m1 + 1 <= 120 then
            List.iter
              (fun (md, da) -> put (Nrmp (nu, m1, m2, da, md)))
              nrmp_args)
        (nrmp_mus nu))
    nrmp_nus;
  List.iter (fun (nu, m1, m2, da, md) -> put (Nrmp (nu, m1, m2, da, md))) nrmp_bad;
  List.rev !out

(* ---- Running the OCaml side ---- *)

let out = open_out_bin "ml.out"
let pr fmt = Printf.fprintf out fmt
let kind_of_id = function 1 -> Pneg | 2 -> Q | 3 -> Ppos | _ -> Pnorm
let mode_of_int = function 1 -> By_x | _ -> By_theta

let emit tag f =
  match f () with
  | a -> pr "%s %s %d\n" tag (hex a.x) a.ix
  | exception Invalid_argument _ -> pr "%s err\n" tag

let vec tag withsig f =
  match f () with
  | v, s ->
    if withsig then pr "%s %d %d\n" tag (Array.length v) s
    else pr "%s %d\n" tag (Array.length v);
    Array.iteri (fun i a -> pr ". %d %s %d\n" (i + 1) (hex a.x) a.ix) v
  | exception Invalid_argument _ -> pr "%s err\n" tag

let run env c =
  match c with
  | Adj (x, ix) -> emit "adj" (fun () -> dxadj env { x; ix })
  | Rt (x, ix) -> emit "rt" (fun () -> dxred env (dxadj env { x; ix }))
  | Add (x, ix, y, iy) ->
    emit "add" (fun () -> dxadd env { x; ix } { x = y; ix = iy })
  | Red (x, ix) -> emit "red" (fun () -> dxred env { x; ix })
  | Con (x, ix) -> emit "con" (fun () -> dxcon env { x; ix })
  | Legf (d, nd, m1, m2, th, id) ->
    vec "legf" false (fun () -> (dxlegf env d nd m1 m2 th (kind_of_id id), 0))
  | Nrmp (nu, m1, m2, da, md) ->
    vec "nrmp" true (fun () -> dxnrmp env nu m1 m2 da (mode_of_int md))

let () =
  List.iteri
    (fun n cfg ->
      let irad, nradpl, dzero, nbits = cfg in
      let cs = cases cfg in
      let ch = open_out_bin (Printf.sprintf "in.%02d" n) in
      Printf.fprintf ch "%d %d %s %d\n" irad nradpl (hex dzero) nbits;
      List.iter (fun c -> Printf.fprintf ch "%s\n" (line c)) cs;
      close_out ch;
      match dxset irad nradpl dzero nbits with
      | env ->
        pr "set ok\n";
        List.iter (run env) cs
      | exception Invalid_argument _ -> pr "set err\n")
    configs;
  close_out out
