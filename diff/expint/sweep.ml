let hex v = Printf.sprintf "%016LX" (Int64.bits_of_float v)

let inp = open_out "sweep.in"
let out = open_out "ml.out"

let emit name x n = Printf.fprintf inp "%s %s %d\n" name (hex x) n

let both x = [ x; -.x ]
let flat l = List.concat_map both l

let small =
  [ 1e-300; 1e-20; 1e-17; 1e-16; 1.1e-16; 1e-15; 1e-12; 1e-9; 1e-8; 1e-6;
    1e-4; 1e-3; 0.01; 0.1 ]

let e1_pts =
  flat
    (small
    @ [ 0.2; 0.3; 0.5; 0.9; 0.999; 1.0; 1.0001; 1.5; 2.0; 3.0; 3.9999; 4.0;
        4.0001; 5.0; 8.0; 10.0; 20.0; 32.0; 32.0001; 50.0; 100.0; 500.0;
        700.0; 701.0; 701.9; 702.0; 710.0; 1000.0 ])

let dli_pts =
  [ 1e-300; 1e-100; 1e-10; 0.001; 0.1; 0.5; 0.9; 0.99; 1.01; 1.1; 1.5; 2.0;
    2.718281828459045; 5.0; 10.0; 100.0; 1e5; 1e10; 1e100; 1e300 ]

let spenc_pts =
  flat
    (small
    @ [ 0.3; 0.5; 0.5000001; 0.7; 0.9; 0.99; 1.0; 1.0001; 1.1; 1.5; 2.0;
        2.0001; 3.0; 5.0; 10.0; 100.0; 1e8; 1e15; 1e16; 1e17; 1e100; 1e300 ])

let daws_pts =
  flat
    (small
    @ [ 0.2; 0.5; 0.9; 1.0; 1.0001; 2.0; 3.0; 4.0; 4.0001; 5.0; 10.0; 1e7;
        1e8; 1e9; 1e10; 1e100; 1e300; 1e307; 1e308 ])

let exprl_pts =
  flat
    (small
    @ [ 0.2; 0.3; 0.4; 0.49; 0.5; 0.500001; 0.6; 1.0; 2.0; 10.0; 100.0; 700.0 ])

let lnrel_pts =
  [ -0.99999; -0.99; -0.9; -0.5; -0.375; -0.3; -0.1; -0.01; -1e-3; -1e-8;
    -1e-16; -1e-300; 0.0; 1e-300; 1e-16; 1e-8; 1e-3; 0.01; 0.1; 0.3; 0.375;
    0.3751; 0.5; 1.0; 10.0; 1e10; 1e300 ]

let cbrt_pts =
  0.0 :: -0.0
  :: flat
       [ 5e-324; 1e-320; 1e-310; min_float; 1e-300; 1e-100; 1e-10; 0.001; 0.125;
         0.5; 1.0; 2.0; 8.0; 27.0; 64.0; 1000.0; 1e10; 1e100; 1e300; max_float ]

let dg_pts =
  0.0
  :: flat
       [ 1e-300; 1e-8; 0.001; 1.0; 1.5; 30.0; 45.0; 60.0; 89.9999; 90.0;
         90.0000001; 135.0; 180.0; 270.0; 360.0; 450.0; 540.0; 1e5; 1e8; 3600.0;
         36000.0; 9e14; 1e15 ]

let atn1_pts =
  flat
    [ 1e-300; 1e-10; 3e-9; 3.4e-9; 1e-8; 1e-5; 0.001; 0.1; 0.5; 0.9; 1.0;
      1.0001; 2.0; 10.0; 1e3; 1e8; 1e15; 1e16; 1.4e16 ]

let ln2r_pts =
  [ -0.99; -0.9; -0.7; -0.625; -0.6249; -0.5; -0.3; -0.1; -1e-3; -1e-8;
    -1e-300; 0.0; 1e-300; 1e-8; 1e-3; 0.1; 0.5; 0.8125; 0.8126; 1.0; 2.0; 10.0;
    1e4; 1e8; 6e8 ]

let pak_ys =
  [ 0.5; 0.75; 1.0; 1.5; 3.0; -0.75; -3.0; 0.0; -0.0; 1e-300; 1e300; 1e-320;
    min_float; max_float ]

let pak_ns = [ 0; 1; -1; 2; -2; 10; -10; 53; -53; 100; -100; 300; -300; 440; -440 ]

let upak_pts =
  0.0 :: -0.0
  :: flat
       [ 5e-324; 1e-320; 1e-310; min_float; 1e-300; 0.25; 0.5; 0.75; 1.0; 1.5;
         2.0; 3.0; 1e10; 1e300; max_float; 0.1; 0.3 ]


let exint_xs =
  [ 0.0; 0.1; 0.3; 0.5; 1.0; 1.5; 1.9999; 2.0; 2.0001; 2.5; 3.0; 4.0; 5.0;
    7.0; 10.0; 20.0; 50.0; 100.0; 500.0; 700.0; 705.0 ]

let exint_ns = [ 1; 2; 3; 4; 5; 10; 25; 99; 101; 200 ]

let () =
  List.iter (fun x -> emit "de1" x 0) e1_pts;
  List.iter (fun x -> emit "dei" x 0) e1_pts;
  List.iter (fun x -> emit "dli" x 0) dli_pts;
  List.iter (fun x -> emit "dspenc" x 0) spenc_pts;
  List.iter (fun x -> emit "ddaws" x 0) daws_pts;
  List.iter (fun x -> emit "dexprl" x 0) exprl_pts;
  List.iter (fun x -> emit "dlnrel" x 0) lnrel_pts;
  List.iter (fun x -> emit "dcbrt" x 0) cbrt_pts;
  List.iter (fun x -> emit "dsindg" x 0) dg_pts;
  List.iter (fun x -> emit "dcosdg" x 0) dg_pts;
  List.iter (fun x -> emit "d9atn1" x 0) atn1_pts;
  List.iter (fun x -> emit "d9ln2r" x 0) ln2r_pts;
  List.iter
    (fun y ->
      let _, ny = Float.frexp y in
      List.iter
        (fun n ->
          let nsum = n + ny in
          if nsum >= -440 && nsum <= 440 then emit "d9pak" y n)
        pak_ns)
    pak_ys;
  List.iter (fun x -> emit "d9upak" x 0) upak_pts;
  List.iter
    (fun x ->
      List.iter
        (fun n ->
          if not (x = 0.0 && n = 1) then begin
            emit "dexint" x n;
            emit "dexint2" x n
          end)
        exint_ns)
    exint_xs;
  close_out inp

open Expint
open Expint.Raw

let p name v = Printf.fprintf out "%s %s\n" name (hex v)

let run_exint tag kode x n =
  let en = Array.make 4 (-1.0) in
  match dexint x n kode 4 1.0e-14 en with
  | nz ->
      Printf.fprintf out "%s_flags %d %d\n" tag nz 0;
      Array.iteri (fun i v -> Printf.fprintf out "%s_%d %s\n" tag (i + 1) (hex v)) en
  | exception Invalid_argument _ -> Printf.fprintf out "%s_flags %d %d\n" tag 0 2

let () =
  List.iter (fun x -> p "de1" (de1 x)) e1_pts;
  List.iter (fun x -> p "dei" (dei x)) e1_pts;
  List.iter (fun x -> p "dli" (dli x)) dli_pts;
  List.iter (fun x -> p "dspenc" (dspenc x)) spenc_pts;
  List.iter (fun x -> p "ddaws" (ddaws x)) daws_pts;
  List.iter (fun x -> p "dexprl" (dexprl x)) exprl_pts;
  List.iter (fun x -> p "dlnrel" (dlnrel x)) lnrel_pts;
  List.iter (fun x -> p "dcbrt" (dcbrt x)) cbrt_pts;
  List.iter (fun x -> p "dsindg" (dsindg x)) dg_pts;
  List.iter (fun x -> p "dcosdg" (dcosdg x)) dg_pts;
  List.iter (fun x -> p "d9atn1" (d9atn1 x)) atn1_pts;
  List.iter (fun x -> p "d9ln2r" (d9ln2r x)) ln2r_pts;
  List.iter
    (fun y ->
      let _, ny = Float.frexp y in
      List.iter
        (fun n ->
          let nsum = n + ny in
          if nsum >= -440 && nsum <= 440 then p "d9pak" (d9pak y n))
        pak_ns)
    pak_ys;
  List.iter
    (fun x ->
      let y, n = d9upak x in
      Printf.fprintf out "d9upak_n %s %d\n" (hex y) n;
      p "d9upak" y)
    upak_pts;
  List.iter
    (fun x ->
      List.iter
        (fun n ->
          if not (x = 0.0 && n = 1) then begin
            run_exint "dexint" 1 x n;
            run_exint "dexint2" 2 x n
          end)
        exint_ns)
    exint_xs;
  close_out out
