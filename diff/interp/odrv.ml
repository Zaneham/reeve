(* Reeve differential, interpolation. Reads cases.hex, which drv.f90 writes,
   and runs the same inputs through Reeve.Interp, printing the results, the
   work arrays and the error count in the same format drv.f90 uses. *)

open Reeve.Interp
open Reeve.Interp.Raw

let nw = 136

let bits s = Int64.of_string ("0x" ^ s)

let dmp tag id idx v =
  Printf.printf "%-8s %6d %6d %016LX\n" tag id idx (Int64.bits_of_float v)

let dmpi tag id idx v =
  Printf.printf "%-8s %6d %6d %016LX\n" tag id idx (Int64.of_int v)

let fill a =
  for i = 1 to nw do
    a.(i - 1) <- -.(float_of_int (1024 + i) /. 8.0)
  done

let read_lines path =
  let ic = open_in path in
  let acc = ref [] in
  (try
     while true do
       acc := input_line ic :: !acc
     done
   with End_of_file -> ());
  close_in ic;
  Array.of_list (List.rev !acc)

let () =
  set_binary_mode_out stdout true;
  let lines = read_lines "cases.hex" in
  let p = ref 0 in
  let next () =
    let s = lines.(!p) in
    incr p;
    s
  in
  let nextf () = Int64.float_of_bits (bits (String.trim (next ()))) in
  let x = Array.make nw 0.0 in
  let y = Array.make nw 0.0 in
  let c = Array.make nw 0.0 in
  let work = Array.make nw 0.0 in
  let yp = Array.make nw 0.0 in
  let d = Array.make nw 0.0 in
  let z = Array.make 3 0.0 in
  while !p < Array.length lines do
    let hdr =
      List.filter (fun s -> s <> "") (String.split_on_char ' ' (next ()))
    in
    let fld k = int_of_string (List.nth hdr k) in
    let id = fld 0 and swp = fld 1 and n = fld 2 and nder = fld 3 in
    let xx = nextf () in
    for i = 0 to 2 do
      z.(i) <- nextf ()
    done;
    Array.fill x 0 nw 0.0;
    Array.fill y 0 nw 0.0;
    for i = 1 to max n 0 do
      x.(i - 1) <- nextf ()
    done;
    for i = 1 to max n 0 do
      y.(i - 1) <- nextf ()
    done;
    let plint () =
      fill c;
      try
        dplint n x y c;
        0
      with Invalid_argument _ -> 1
    in
    match swp with
    | 1 ->
      let xerr = plint () in
      for i = 1 to max n 1 do
        dmp "a.c" id i c.(i - 1)
      done;
      dmpi "a.xerr" id 0 xerr
    | 2 ->
      ignore (plint ());
      let nd1 = max nder 1 in
      fill work;
      fill yp;
      let yfit, ierr = dpolvl nder xx yp n x c work in
      dmp "b.yfit" id 0 yfit;
      dmpi "b.ierr" id 0 ierr;
      for i = 1 to nd1 do
        dmp "b.yp" id i yp.(i - 1)
      done;
      for i = 1 to 2 * n do
        dmp "b.work" id i work.(i - 1)
      done
    | 3 ->
      ignore (plint ());
      fill work;
      fill d;
      dpolcf xx n x c d work;
      for i = 1 to n do
        dmp "c.d" id i d.(i - 1)
      done;
      for i = 1 to 2 * n do
        dmp "c.work" id i work.(i - 1)
      done;
      for k = 1 to 3 do
        let acc = ref d.(n - 1) in
        for i = n - 1 downto 1 do
          acc := (!acc *. (z.(k - 1) -. xx)) +. d.(i - 1)
        done;
        dmp "c.horn" id k !acc;
        fill work;
        fill yp;
        let yfit, _ = dpolvl 0 z.(k - 1) yp n x c work in
        dmp "c.yfit" id k yfit
      done
    | _ -> failwith "odrv: unknown sweep"
  done
