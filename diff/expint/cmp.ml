let hex v = Printf.sprintf "%016LX" (Int64.bits_of_float v)
let unhex s = Int64.float_of_bits (Int64.of_string ("0x" ^ s))

let pts =
  let l = ref [] in
  let add x = l := x :: -.x :: !l in
  for k = 0 to 80 do
    add (Float.ldexp 1.0 (-k));
    add (Float.ldexp 1.3 (-k));
    add (Float.ldexp 1.9999 (-k))
  done;
  List.iter add
    [ 0.375; 0.3750001; 0.3749999; 0.5; 0.5000001; 0.4999999; 0.8125; 1.0;
      2.0; 8.0; 27.0; 1e-17; 1e-16; 1e-15; 1e-12; 1e-9; 1e-6; 1e-3; 0.01; 0.1;
      0.2; 0.3; 0.49; 0.7; 0.9; 1.5; 3.0; 10.0; 100.0 ];
  let st = Random.State.make [| 7 |] in
  for _ = 1 to 4000 do
    let e = Random.State.int st 70 in
    let m = 1.0 +. Random.State.float st 1.0 in
    add (Float.ldexp m (-e))
  done;
  !l

let () =
  match Sys.argv.(1) with
  | "gen" ->
      let oc = open_out "qref.in" in
      List.iter
        (fun x ->
          if x > -1.0 then Printf.fprintf oc "lnrel %s\n" (hex x);
          Printf.fprintf oc "exprl %s\n" (hex x);
          Printf.fprintf oc "cbrt %s\n" (hex x))
        pts;
      close_out oc
  | "cmp" ->
      let ic = open_in "qref.out" in
      let tbl = Hashtbl.create 100 in
      (try
         while true do
           let line = input_line ic in
           match String.split_on_char ' ' line with
           | [ nm; xh; hih; loh ] ->
               let cur = try Hashtbl.find tbl nm with Not_found -> [] in
               Hashtbl.replace tbl nm ((unhex xh, unhex hih, unhex loh) :: cur)
           | _ -> ()
         done
       with End_of_file -> ());
      close_in ic;
      let ulps v hi lo =
        if hi = 0.0 then if v = 0.0 then 0.0 else infinity
        else
          let _, e = Float.frexp hi in
          let u = Float.ldexp 1.0 (e - 53) in
          (v -. hi -. lo) /. u
      in
      let report nm f g =
        let rows = try Hashtbl.find tbl nm with Not_found -> [] in
        let worst_f = ref 0.0
        and worst_g = ref 0.0
        and at_f = ref 0.0
        and at_g = ref 0.0
        and n = ref 0
        and fwin = ref 0
        and gwin = ref 0
        and tie = ref 0 in
        List.iter
          (fun (x, hi, lo) ->
            let a = Float.abs (ulps (f x) hi lo) in
            let b = Float.abs (ulps (g x) hi lo) in
            if Float.is_nan a || Float.is_nan b then ()
            else begin
              incr n;
              if a > !worst_f then begin
                worst_f := a;
                at_f := x
              end;
              if b > !worst_g then begin
                worst_g := b;
                at_g := x
              end;
              if a < b then incr fwin
              else if b < a then incr gwin
              else incr tie
            end)
          rows;
        Printf.printf
          "%-6s n=%d  slatec max %.3f ulp (at %.17g)  stdlib max %.3f ulp (at %.17g)  slatec closer %d, stdlib closer %d, equal %d\n"
          nm !n !worst_f !at_f !worst_g !at_g !fwin !gwin !tie
      in
      report "exprl"
        (fun x -> Expint.dexprl x)
        (fun x -> if x = 0.0 then 1.0 else Float.expm1 x /. x);
      report "lnrel" (fun x -> Expint.dlnrel x) (fun x -> Float.log1p x);
      report "cbrt" (fun x -> Expint.dcbrt x) (fun x -> Float.cbrt x);
      (* exact-identity counts for the pack pair *)
      let pak_diff = ref 0 and pak_total = ref 0 and first = ref None in
      List.iter
        (fun y ->
          for n = -1200 to 1200 do
            incr pak_total;
            let a = try Some (Expint.d9pak y n) with Invalid_argument _ -> None in
            let b = Float.ldexp y n in
            let same =
              match a with
              | None -> false
              | Some v -> Int64.bits_of_float v = Int64.bits_of_float b
            in
            if not same then begin
              incr pak_diff;
              if !first = None then first := Some (y, n, a, b)
            end
          done)
        [ 0.5; 0.75; 1.0; 1.5; 3.0; -0.75; 1e-300; 1e300; 0.0 ];
      Printf.printf "d9pak vs ldexp: %d of %d (y,n) pairs differ\n" !pak_diff
        !pak_total;
      (match !first with
      | Some (y, n, a, b) ->
          Printf.printf "  first: d9pak %.17g %d = %s, ldexp = %s\n" y n
            (match a with None -> "raise" | Some v -> hex v)
            (hex b)
      | None -> ());
      (* where exactly do they part company *)
      let y = 0.75 in
      let probe n =
        let a = try Some (Expint.d9pak y n) with Invalid_argument _ -> None in
        let b = Float.ldexp y n in
        Printf.printf "  n=%5d d9pak=%-18s ldexp=%s\n" n
          (match a with None -> "raise" | Some v -> hex v)
          (hex b)
      in
      List.iter probe [ -1023; -1022; -1021; -1020; 1023; 1024; 1025 ];
      let cz = Expint.dcbrt (-0.0) in
      Printf.printf "dcbrt (-0.) = %s, Float.cbrt (-0.) = %s\n" (hex cz)
        (hex (Float.cbrt (-0.0)))
  | _ -> ()
