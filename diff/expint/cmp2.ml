let unhex s = Int64.float_of_bits (Int64.of_string ("0x" ^ s))

let () =
  let ic = open_in "qref.out" in
  let tbl = Hashtbl.create 10 in
  (try
     while true do
       match String.split_on_char ' ' (input_line ic) with
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
      (v -. hi -. lo) /. Float.ldexp 1.0 (e - 53)
  in
  let band nm lo hi variants =
    let rows = List.filter (fun (x, _, _) ->
        let a = Float.abs x in a >= lo && a < hi)
        (try Hashtbl.find tbl nm with Not_found -> []) in
    let n = List.length rows in
    if n > 0 then begin
      Printf.printf "%-6s |x| in [%g, %g)  n=%d\n" nm lo hi n;
      List.iter (fun (label, f) ->
          let w = ref 0.0 and s = ref 0.0 and cnt = ref 0 in
          List.iter (fun (x, h, l) ->
              let e = Float.abs (ulps (f x) h l) in
              if Float.is_finite e then begin
                if e > !w then w := e;
                s := !s +. (e *. e);
                incr cnt
              end) rows;
          Printf.printf "    %-22s max %10.3f ulp   rms %8.3f ulp\n" label !w
            (sqrt (!s /. float_of_int (max 1 !cnt)))) variants
    end
  in
  let bands = [ (0.0, 1e-8); (1e-8, 1e-4); (1e-4, 0.1); (0.1, 1e3) ] in
  List.iter (fun (lo, hi) ->
      band "exprl" lo hi
        [ ("slatec dexprl", Expint.dexprl);
          ("stdlib expm1 x /. x", fun x -> if x = 0.0 then 1.0 else Float.expm1 x /. x);
          ("naive (exp x-1)/x", fun x -> if x = 0.0 then 1.0 else (exp x -. 1.0) /. x) ])
    bands;
  List.iter (fun (lo, hi) ->
      band "lnrel" lo hi
        [ ("slatec dlnrel", Expint.dlnrel);
          ("stdlib log1p", Float.log1p);
          ("naive log (1+x)", fun x -> log (1.0 +. x)) ])
    bands;
  List.iter (fun (lo, hi) ->
      band "cbrt" lo hi
        [ ("slatec dcbrt", Expint.dcbrt); ("stdlib cbrt", Float.cbrt) ])
    bands
