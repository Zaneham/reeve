open Reeve.Sort
open Reeve.Sort.Raw

let gen n sel =
  let dx = Array.make (max n 1) 0.0 in
  let s = ref (123456789 + (sel * 7)) in
  for i = 0 to n - 1 do
    s := (1103515245 * !s + 12345) mod 2147483648;
    dx.(i) <-
      (if sel = 1 then float_of_int ((!s / 65536 mod 1000) - 500) /. 8.0
       else if sel = 2 then float_of_int (!s / 65536 mod 7)
       else if sel = 3 then 2.5
       else if sel = 4 then float_of_int (i + 1) /. 4.0
       else if sel = 5 then float_of_int (n - i) /. 4.0
       else float_of_int (min (i + 1) (n - i)) /. 2.0)
  done;
  dx

let dumpd tag n (a : float array) =
  for i = 1 to n do
    Printf.printf "%s %6d %016LX\n" tag i (Int64.bits_of_float a.(i - 1))
  done

let dumpi tag n (a : int array) =
  for i = 1 to n do
    Printf.printf "%s %6d %8d\n" tag i (a.(i - 1) + 1)
  done

let flag_of k =
  match k with
  | 1 -> Increasing
  | 2 -> Increasing_carry
  | -1 -> Decreasing
  | -2 -> Decreasing_carry
  | _ -> assert false

let nmax = 2000

let () =
  let sizes = [ 1; 2; 3; 5; 8; 11; 12; 13; 100; 257; 1000; 2000 ] in
  let flags = [ 1; 2; -1; -2 ] in
  for sel = 1 to 6 do
    List.iter
      (fun n ->
        List.iter
          (fun kflag ->
            let dx = gen n sel in
            let dy = Array.init n (fun p -> float_of_int (p + 1)) in
            dsort dx dy n (flag_of kflag);
            let tag = Printf.sprintf "dsort s%d n%5d k%2d" sel n kflag in
            dumpd (tag ^ " x") n dx;
            dumpd (tag ^ " y") n dy;
            let dx = gen n sel in
            let iperm = Array.make n 0 in
            dpsort dx n iperm (flag_of kflag);
            let tag = Printf.sprintf "dpsort s%d n%5d k%2d" sel n kflag in
            dumpd (tag ^ " x") n dx;
            dumpi (tag ^ " p") n iperm;
            Printf.printf "%s ier %3d\n" tag 0)
          flags;
        let dx = gen n sel in
        let iperm = Array.make n 0 in
        dpsort dx n iperm Increasing;
        let dx = gen n sel in
        dpperm dx n iperm;
        let tag = Printf.sprintf "dpperm s%d n%5d" sel n in
        dumpd (tag ^ " x") n dx;
        dumpi (tag ^ " p") n iperm;
        Printf.printf "%s ier %3d\n" tag 0;
        let m1 = n / 2 in
        let m2 = n - m1 in
        let tc = Array.make (3 * nmax) (-7.5) in
        let g = gen n sel in
        Array.blit g 0 tc 0 n;
        if m1 > 0 then begin
          let a = Array.sub tc 0 m1 in
          dsort a [||] m1 Increasing;
          Array.blit a 0 tc 0 m1
        end;
        let b = Array.sub tc m1 m2 in
        dsort b [||] m2 Increasing;
        Array.blit b 0 tc m1 m2;
        d1merg tc 0 m1 m1 m2 n;
        let tag = Printf.sprintf "d1merg s%d n%5d" sel n in
        dumpd (tag ^ " t") (2 * n) tc)
      sizes
  done
