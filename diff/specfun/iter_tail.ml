
let ends_with s suf =
  let n = String.length s and m = String.length suf in
  n >= m && String.sub s (n - m) m = suf

let hits = ref 0

let run f =
  try ignore (f () : float) with
  | Invalid_argument m -> if ends_with m "iteration limit reached" then incr hits
  | _ -> ()

let () =
  let vals = [| 0.0; drf_lolim; 1e-300; 1e-10; 1.0; 1e10; 1e300; drf_uplim |] in
  Array.iter
    (fun x ->
      Array.iter
        (fun y ->
          run (fun () -> drc x y);
          Array.iter
            (fun z ->
              run (fun () -> drf x y z);
              run (fun () -> drd x y z))
            vals)
        vals)
    vals;
  Printf.printf "maxit=%d  iteration-limit hits=%d\n" carlson_maxit !hits
