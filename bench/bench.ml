(* Reeve, Copyright 2026 Zane Hambly.
   One executable that times the float-heavy paths, so the same numbers can be
   taken under two compilers and compared. *)

let floor_seconds = 0.5
let repeats = 1
let sink = ref 0.0
let now () = Unix.gettimeofday ()

(* The checksum runs once, outside the timing, over a fixed single pass, so it
   does not move when the calibration picks a different iteration count. *)
let bench name check run =
  let csum = check () in
  let once n = let t = now () in run n; now () -. t in
  let n = ref 1 in
  while once !n < floor_seconds do n := !n * 2 done;
  let best = ref infinity and worst = ref 0.0 in
  for _ = 1 to repeats do
    let t = once !n in
    if t < !best then best := t;
    if t > !worst then worst := t
  done;
  Printf.printf "%-13s %9d %10.2f %7.2f %.17g\n" name !n
    (!best /. float_of_int !n *. 1.0e9)
    ((!worst -. !best) /. !best *. 100.0)
    csum

let sweep lo hi k =
  Array.init k (fun i ->
      lo +. ((hi -. lo) *. float_of_int i /. float_of_int (k - 1)))

let scalar name f lo hi =
  let xs = sweep lo hi 1000 in
  let pass () =
    let acc = ref 0.0 in
    for i = 0 to Array.length xs - 1 do
      acc := !acc +. f (Sys.opaque_identity xs.(i))
    done;
    !acc
  in
  bench name pass (fun n -> for _ = 1 to n do sink := !sink +. pass () done)

let vec k = Array.init k (fun i -> 1.0 +. (float_of_int i *. 1.0e-3))

let spd k =
  let a = Array.make (k * k) 0.0 in
  for j = 0 to k - 1 do
    for i = 0 to k - 1 do
      a.(i + (j * k)) <- 1.0 /. (1.0 +. float_of_int (abs (i - j)))
    done;
    a.(j + (j * k)) <- float_of_int k
  done;
  a

let noise k =
  let g = ref 123456789 in
  Array.init k (fun _ ->
      g := ((!g * 1103515245) + 12345) land 0x3FFFFFFF;
      float_of_int !g)

let () =
  Printf.printf "%-13s %9s %10s %7s %s\n" "bench" "iters" "ns/iter" "spread"
    "checksum";

  scalar "dbesi0" Reeve.Bessel.dbesi0 0.5 20.0;
  scalar "dbesk0" Reeve.Bessel.dbesk0 0.5 20.0;
  scalar "dgamln" Reeve.Specfun.dgamln 0.5 60.0;
  scalar "dpsi" Reeve.Specfun.dpsi 0.5 60.0;
  scalar "dai" Reeve.Specfun.dai (-8.0) 4.0;
  scalar "de1" Reeve.Expint.de1 0.25 30.0;
  scalar "dexprl" Reeve.Expint.dexprl (-2.0) 2.0;
  scalar "dlnrel" Reeve.Expint.dlnrel (-0.5) 4.0;
  scalar "dgamr" Reeve.Gamma.dgamr 0.5 60.0;
  scalar "dcot" Reeve.Specfun.Raw.dcot 0.5 60.0;
  scalar "dsindg" Reeve.Expint.dsindg (-720.0) 720.0;
  scalar "dcosdg" Reeve.Expint.dcosdg (-720.0) 720.0;

  let xs = sweep 0.25 8.0 1000 in
  let rf () =
    let acc = ref 0.0 in
    for i = 0 to Array.length xs - 1 do
      acc := !acc +. Reeve.Specfun.drf xs.(i) (xs.(i) +. 1.0) (xs.(i) +. 2.0)
    done;
    !acc
  in
  bench "drf" rf (fun n -> for _ = 1 to n do sink := !sink +. rf () done);

  let jx = sweep 0.5 40.0 200 and jy = Array.make 8 0.0 in
  let besj () =
    let acc = ref 0.0 in
    for i = 0 to Array.length jx - 1 do
      ignore (Reeve.Bessel.Raw.dbesj jx.(i) 0.0 8 jy);
      acc := !acc +. jy.(7)
    done;
    !acc
  in
  bench "dbesj.seq8" besj
    (fun n -> for _ = 1 to n do sink := !sink +. besj () done);

  let ax = vec 10_000 and ay = vec 10_000 in
  bench "daxpy.10k"
    (fun () ->
      let y = Array.copy ay in
      Reeve.Blas.Raw.daxpy 10_000 1.000001 ax 0 1 y 0 1;
      y.(0) +. y.(9_999))
    (fun n ->
      for _ = 1 to n do
        Reeve.Blas.Raw.daxpy 10_000 1.000001 ax 0 1 ay 0 1
      done;
      sink := !sink +. ay.(0));

  let dx = vec 10_000 and dy = vec 10_000 in
  bench "ddot.10k"
    (fun () -> Reeve.Blas.Raw.ddot 10_000 dx 0 1 dy 0 1)
    (fun n ->
      for _ = 1 to n do
        sink := !sink +. Reeve.Blas.Raw.ddot 10_000 dx 0 1 dy 0 1
      done);

  bench "dnrm2.10k"
    (fun () -> Reeve.Blas.Raw.dnrm2 10_000 dx 0 1)
    (fun n ->
      for _ = 1 to n do
        sink := !sink +. Reeve.Blas.Raw.dnrm2 10_000 dx 0 1
      done);

  let k = 128 in
  let ga = spd k and gb = spd k and gc = Array.make (k * k) 0.0 in
  let gemm () =
    Reeve.Blasmat.Raw.dgemm Reeve.Blasmat.No_trans Reeve.Blasmat.No_trans k k k
      1.0 ga 0 k gb 0 k 0.0 gc 0 k;
    let s = ref 0.0 in
    Array.iter (fun v -> s := !s +. v) gc;
    !s
  in
  bench "dgemm.128" gemm
    (fun n ->
      for _ = 1 to n do
        Reeve.Blasmat.Raw.dgemm Reeve.Blasmat.No_trans Reeve.Blasmat.No_trans k
          k k 1.0 ga 0 k gb 0 k 0.0 gc 0 k
      done;
      sink := !sink +. gc.(0));

  let pristine = spd k and work = Array.make (k * k) 0.0 in
  let ipiv = Array.make k 0 in
  bench "dgetrf.128"
    (fun () ->
      Array.blit pristine 0 work 0 (k * k);
      ignore (Reeve.Lapack.Raw.dgetrf k k work 0 k ipiv 0);
      let s = ref 0.0 in
      Array.iter (fun v -> s := !s +. v) work;
      !s)
    (fun n ->
      for _ = 1 to n do
        Array.blit pristine 0 work 0 (k * k);
        ignore (Reeve.Lapack.Raw.dgetrf k k work 0 k ipiv 0)
      done;
      sink := !sink +. work.(0));

  bench "dpotrf.128"
    (fun () ->
      Array.blit pristine 0 work 0 (k * k);
      ignore (Reeve.Lapack.Raw.dpotrf Reeve.Blasmat.Upper k work 0 k);
      let s = ref 0.0 in
      Array.iter (fun v -> s := !s +. v) work;
      !s)
    (fun n ->
      for _ = 1 to n do
        Array.blit pristine 0 work 0 (k * k);
        ignore (Reeve.Lapack.Raw.dpotrf Reeve.Blasmat.Upper k work 0 k)
      done;
      sink := !sink +. work.(0));

  let raw = noise 10_000 and sorted = Array.make 10_000 0.0 in
  bench "dsort.10k"
    (fun () ->
      Array.blit raw 0 sorted 0 10_000;
      Reeve.Sort.Raw.dsort sorted [||] 10_000 Reeve.Sort.Raw.Increasing;
      sorted.(0) +. sorted.(9_999))
    (fun n ->
      for _ = 1 to n do
        Array.blit raw 0 sorted 0 10_000;
        Reeve.Sort.Raw.dsort sorted [||] 10_000 Reeve.Sort.Raw.Increasing
      done;
      sink := !sink +. sorted.(0));

  if Float.is_nan !sink then print_string "sink went NaN\n"
