(* SLATEC in OCaml, Copyright 2026 Zane Hambly.
   Checks for the data handling routines. A sort has no reference number to
   check against, so each case checks the result against the problem: the
   output is ordered, it is a permutation of the input, the carried array
   moved with it, and the permutation routines agree with the permutation the
   sort reported. *)

open Reeve.Sort.Raw
module Sort = Reeve.Sort

let bad = ref 0

let fail what =
  Printf.printf "FAIL %s\n" what;
  bad := !bad + 1

let check what c = if not c then fail what

let raises what f =
  match f () with
  | () -> fail (what ^ " did not raise")
  | exception Invalid_argument _ -> ()

let ordered what up (a : float array) n =
  for p = 0 to n - 2 do
    let ok = if up then a.(p) <= a.(p + 1) else a.(p) >= a.(p + 1) in
    if not ok then fail (Printf.sprintf "%s order at %d" what p)
  done

let is_perm what (p : int array) n =
  let seen = Array.make (max n 1) false in
  for q = 0 to n - 1 do
    let v = p.(q) in
    if v < 0 || v >= n then fail (Printf.sprintf "%s range at %d" what q)
    else if seen.(v) then fail (Printf.sprintf "%s repeat at %d" what q)
    else seen.(v) <- true
  done

let count x (a : float array) n =
  let c = ref 0 in
  for p = 0 to n - 1 do
    if a.(p) = x then c := !c + 1
  done;
  !c

let multiset what (a : float array) b n =
  for p = 0 to n - 1 do
    if count a.(p) a n <> count a.(p) b n then
      fail (Printf.sprintf "%s multiset at %d" what p)
  done

let exact what (a : float array) (b : float array) n =
  for p = 0 to n - 1 do
    if a.(p) <> b.(p) then fail (Printf.sprintf "%s[%d]" what p)
  done

let exact_ints what (a : int array) (b : int array) n =
  for p = 0 to n - 1 do
    if a.(p) <> b.(p) then fail (Printf.sprintf "%s[%d]" what p)
  done

let tags n = Array.init n (fun p -> float_of_int p)

let paired what (orig : float array) (a : float array) (tag : float array) n =
  let idx = Array.init n (fun p -> int_of_float tag.(p)) in
  is_perm (what ^ " tag") idx n;
  for p = 0 to n - 1 do
    if a.(p) <> orig.(idx.(p)) then
      fail (Printf.sprintf "%s pair at %d" what p)
  done

let picked what (orig : float array) (a : float array) (ip : int array) n =
  for p = 0 to n - 1 do
    if a.(p) <> orig.(ip.(p)) then
      fail (Printf.sprintf "%s pick at %d" what p)
  done

(* ---- dsort ---- *)

let dsort_case what xs =
  let n = Array.length xs in
  let orig = Array.copy xs in
  let a = Array.copy xs in
  let y = tags n in
  dsort a y n Increasing_carry;
  ordered (what ^ " up carry") true a n;
  paired (what ^ " up carry") orig a y n;
  let a = Array.copy xs in
  let y = tags n in
  dsort a y n Decreasing_carry;
  ordered (what ^ " down carry") false a n;
  paired (what ^ " down carry") orig a y n;
  let a = Array.copy xs in
  let y = tags n in
  dsort a y n Increasing;
  ordered (what ^ " up") true a n;
  multiset (what ^ " up") orig a n;
  exact (what ^ " up dy held") y (tags n) n;
  let a = Array.copy xs in
  let y = tags n in
  dsort a y n Decreasing;
  ordered (what ^ " down") false a n;
  multiset (what ^ " down") orig a n;
  exact (what ^ " down dy held") y (tags n) n

(* ---- dpsort and dpperm ---- *)

let dpsort_case what xs =
  let n = Array.length xs in
  let orig = Array.copy xs in
  let a = Array.copy xs in
  let ip = Array.make n 0 in
  dpsort a n ip Increasing;
  is_perm (what ^ " up perm") ip n;
  exact (what ^ " up dx held") a orig n;
  ordered (what ^ " up order") true (Array.init n (fun p -> orig.(ip.(p)))) n;
  let keep = Array.copy ip in
  let a = Array.copy xs in
  let ip = Array.make n 0 in
  dpsort a n ip Increasing_carry;
  is_perm (what ^ " up carry perm") ip n;
  ordered (what ^ " up carry order") true a n;
  picked (what ^ " up carry") orig a ip n;
  let a = Array.copy xs in
  let ip = Array.make n 0 in
  dpsort a n ip Decreasing;
  is_perm (what ^ " down perm") ip n;
  exact (what ^ " down dx held") a orig n;
  ordered (what ^ " down order") false (Array.init n (fun p -> orig.(ip.(p)))) n;
  let a = Array.copy xs in
  let ip = Array.make n 0 in
  dpsort a n ip Decreasing_carry;
  is_perm (what ^ " down carry perm") ip n;
  ordered (what ^ " down carry order") false a n;
  picked (what ^ " down carry") orig a ip n;
  let b = Array.copy xs in
  let ip = Array.copy keep in
  dpperm b n ip;
  exact_ints (what ^ " dpperm iperm held") ip keep n;
  ordered (what ^ " dpperm order") true b n;
  picked (what ^ " dpperm") orig b keep n;
  let t = tags n in
  dpperm t n ip;
  exact_ints
    (what ^ " dpperm carried")
    (Array.init n (fun p -> int_of_float t.(p)))
    keep n

let both_cases what xs =
  dsort_case what xs;
  dpsort_case what xs

(* ---- Inputs ---- *)

let one () = both_cases "one" [| 3.5 |]
let two_up () = both_cases "two up" [| 1.0; 2.0 |]
let two_down () = both_cases "two down" [| 2.0; 1.0 |]
let two_equal () = both_cases "two equal" [| 2.0; 2.0 |]

let all_equal () =
  both_cases "all equal" (Array.make 17 (-4.25))

let already_sorted () =
  both_cases "already sorted" (Array.init 64 (fun p -> float_of_int p))

let reverse_sorted () =
  both_cases "reverse sorted" (Array.init 64 (fun p -> float_of_int (63 - p)))

let saw_tooth () =
  both_cases "saw tooth"
    (Array.init 51 (fun p -> float_of_int (p mod 7) -. 3.0))

let two_values () =
  both_cases "two values"
    (Array.init 40 (fun p -> if p mod 2 = 0 then 1.0 else 0.0))

let organ_pipe () =
  both_cases "organ pipe"
    (Array.init 65 (fun p -> float_of_int (min p (64 - p))))

let signed_mix () =
  both_cases "signed mix"
    [| 0.0; -0.0; 1.0; -1.0; infinity; neg_infinity; 1e-300; -1e-300; 5.0 |]

let random_sizes () =
  let st = Random.State.make [| 20261001 |] in
  List.iter
    (fun n ->
      let xs = Array.init n (fun _ -> Random.State.float st 2.0 -. 1.0) in
      both_cases (Printf.sprintf "random %d" n) xs;
      let ys =
        Array.init n (fun _ -> float_of_int (Random.State.int st 5))
      in
      both_cases (Printf.sprintf "random ties %d" n) ys)
    [ 2; 3; 4; 5; 7; 11; 12; 13; 21; 22; 33; 100; 257; 1000 ]

let tail_untouched () =
  let n = 9 in
  let a = Array.init 16 (fun p -> float_of_int (16 - p)) in
  let y = Array.init 16 (fun p -> float_of_int p) in
  let orig = Array.copy a in
  dsort a y n Increasing_carry;
  ordered "tail dsort order" true a n;
  for p = n to 15 do
    if a.(p) <> orig.(p) then fail (Printf.sprintf "tail dsort dx[%d]" p);
    if y.(p) <> float_of_int p then
      fail (Printf.sprintf "tail dsort dy[%d]" p)
  done;
  let b = Array.copy orig in
  let ip = Array.make 16 (-1) in
  dpsort b n ip Increasing_carry;
  ordered "tail dpsort order" true b n;
  is_perm "tail dpsort perm" ip n;
  for p = n to 15 do
    if b.(p) <> orig.(p) then fail (Printf.sprintf "tail dpsort dx[%d]" p);
    if ip.(p) <> -1 then fail (Printf.sprintf "tail dpsort iperm[%d]" p)
  done;
  let c = Array.copy orig in
  dpperm c n ip;
  ordered "tail dpperm order" true c n;
  for p = n to 15 do
    if c.(p) <> orig.(p) then fail (Printf.sprintf "tail dpperm dx[%d]" p)
  done

let dsort_ignores_dy () =
  let a = [| 3.0; 1.0; 2.0 |] in
  dsort a [||] 3 Increasing;
  ordered "empty dy" true a 3;
  let b = [| 3.0; 1.0; 2.0 |] in
  dsort b [||] 3 Decreasing;
  ordered "empty dy down" false b 3

let dpperm_identity () =
  let a = [| 4.0; 1.0; 9.0; 2.0 |] in
  let ip = [| 0; 1; 2; 3 |] in
  dpperm a 4 ip;
  exact "dpperm identity" a [| 4.0; 1.0; 9.0; 2.0 |] 4;
  exact_ints "dpperm identity iperm" ip [| 0; 1; 2; 3 |] 4

let dpperm_cycle () =
  let a = [| 10.0; 20.0; 30.0; 40.0; 50.0 |] in
  let ip = [| 4; 0; 1; 2; 3 |] in
  dpperm a 5 ip;
  exact "dpperm cycle" a [| 50.0; 10.0; 20.0; 30.0; 40.0 |] 5;
  exact_ints "dpperm cycle iperm" ip [| 4; 0; 1; 2; 3 |] 5

let dpperm_swaps () =
  let a = [| 1.0; 2.0; 3.0; 4.0 |] in
  let ip = [| 1; 0; 3; 2 |] in
  dpperm a 4 ip;
  exact "dpperm swaps" a [| 2.0; 1.0; 4.0; 3.0 |] 4

let errors () =
  raises "dsort n zero" (fun () -> dsort [| 1.0 |] [| 1.0 |] 0 Increasing);
  raises "dsort n negative" (fun () -> dsort [| 1.0 |] [| 1.0 |] (-1) Increasing);
  raises "dpsort n zero" (fun () -> dpsort [| 1.0 |] 0 [| 0 |] Increasing);
  raises "dpperm n zero" (fun () -> dpperm [| 1.0 |] 0 [| 0 |]);
  raises "dpperm repeat" (fun () -> dpperm [| 1.0; 2.0 |] 2 [| 0; 0 |]);
  raises "dpperm high" (fun () -> dpperm [| 1.0; 2.0 |] 2 [| 0; 2 |]);
  raises "dpperm low" (fun () -> dpperm [| 1.0; 2.0 |] 2 [| -1; 1 |]);
  raises "dpperm negative" (fun () -> dpperm [| 1.0; 2.0 |] 2 [| 1; -2 |])

(* ---- d1merg ---- *)

let merge_case what xs ys =
  let m1 = Array.length xs in
  let m2 = Array.length ys in
  let n = m1 + m2 in
  let tcos = Array.make (n + n) 0.0 in
  Array.blit xs 0 tcos 0 m1;
  Array.blit ys 0 tcos m1 m2;
  let src = Array.append xs ys in
  d1merg tcos 0 m1 m1 m2 n;
  let got = Array.sub tcos n n in
  ordered (what ^ " order") true got n;
  multiset (what ^ " multiset") src got n;
  exact (what ^ " sources held") (Array.sub tcos 0 n) src n

let merge_plain () =
  merge_case "merge plain" [| 1.0; 4.0; 9.0 |] [| 2.0; 3.0; 10.0; 11.0 |]

let merge_interleaved () =
  merge_case "merge interleaved"
    (Array.init 8 (fun p -> float_of_int (2 * p)))
    (Array.init 8 (fun p -> float_of_int ((2 * p) + 1)))

let merge_disjoint_ranges () =
  merge_case "merge disjoint ranges" [| 1.0; 2.0; 3.0 |] [| 7.0; 8.0 |];
  merge_case "merge disjoint ranges back" [| 7.0; 8.0 |] [| 1.0; 2.0; 3.0 |]

let merge_ties () =
  merge_case "merge ties" [| 1.0; 2.0; 2.0; 5.0 |] [| 2.0; 2.0; 5.0 |];
  merge_case "merge all equal" (Array.make 5 3.0) (Array.make 4 3.0)

let merge_singletons () =
  merge_case "merge one one" [| 2.0 |] [| 1.0 |];
  merge_case "merge one many" [| 6.0 |] [| 1.0; 2.0; 7.0 |]

let merge_empty () =
  let tcos = [| 1.0; 2.0; 3.0; 0.0; 0.0; 0.0 |] in
  d1merg tcos 0 3 3 0 3;
  exact "merge m2 zero" tcos [| 1.0; 2.0; 3.0; 1.0; 2.0; 3.0 |] 6;
  let tcos = [| 0.0; 0.0; 4.0; 5.0; 0.0; 0.0 |] in
  d1merg tcos 0 0 2 2 4;
  exact "merge m1 zero" tcos [| 0.0; 0.0; 4.0; 5.0; 4.0; 5.0 |] 6;
  let tcos = [| 1.0; 2.0; 3.0 |] in
  d1merg tcos 0 0 1 0 2;
  exact "merge both zero" tcos [| 1.0; 2.0; 3.0 |] 3

let merge_offsets () =
  let tcos = Array.make 16 (-1.0) in
  tcos.(2) <- 1.0;
  tcos.(3) <- 5.0;
  tcos.(4) <- 6.0;
  tcos.(8) <- 2.0;
  tcos.(9) <- 4.0;
  d1merg tcos 2 3 8 2 10;
  exact "merge offsets"
    (Array.sub tcos 10 5)
    [| 1.0; 2.0; 4.0; 5.0; 6.0 |]
    5;
  check "merge offsets tail" (tcos.(15) = -1.0);
  check "merge offsets gap" (tcos.(5) = -1.0 && tcos.(7) = -1.0)

let merge_random () =
  let st = Random.State.make [| 20261002 |] in
  let run n =
    let a = Array.make n 0.0 in
    let v = ref 0.0 in
    for p = 0 to n - 1 do
      v := !v +. float_of_int (Random.State.int st 3);
      a.(p) <- !v
    done;
    a
  in
  List.iter
    (fun (m1, m2) ->
      merge_case (Printf.sprintf "merge random %d %d" m1 m2) (run m1) (run m2))
    [ (1, 1); (1, 9); (9, 1); (5, 5); (17, 3); (3, 17); (40, 40) ]

(* ---- The OCaml surface ---- *)

let surface_sorts () =
  let a = [| 3.0; -1.0; 2.0; 0.0 |] in
  Sort.sort a;
  exact "surface sort up" a [| -1.0; 0.0; 2.0; 3.0 |] 4;
  Sort.sort ~desc:true a;
  exact "surface sort down" a [| 3.0; 2.0; 0.0; -1.0 |] 4

let surface_carries () =
  let a = [| 3.0; 1.0; 2.0 |] and y = [| 30.0; 10.0; 20.0 |] in
  Sort.sort ~carry:y a;
  exact "surface carry keys" a [| 1.0; 2.0; 3.0 |] 3;
  exact "surface carry values" y [| 10.0; 20.0; 30.0 |] 3;
  let b = [| 1.0; 2.0 |] and z = [| 10.0; 20.0; 99.0 |] in
  Sort.sort ~desc:true ~carry:z b;
  exact "surface carry down keys" b [| 2.0; 1.0 |] 2;
  exact "surface carry down values" z [| 20.0; 10.0; 99.0 |] 3

let surface_indexes () =
  let a = [| 5.0; 2.0; 9.0; 4.0; 2.0 |] in
  let keep = Array.copy a in
  let p = Sort.sort_index a in
  exact "surface index leaves the keys" a keep 5;
  is_perm "surface index is a permutation" p 5;
  let got = Array.map (fun i -> a.(i)) p in
  ordered "surface index order" true got 5;
  multiset "surface index multiset" got a 5;
  let b = Array.copy a in
  Sort.permute b p;
  exact "surface permute replays the index" b got 5;
  let q = Sort.sort_index ~desc:true a in
  exact "surface index down leaves the keys" a keep 5;
  ordered "surface index down order" false (Array.map (fun i -> a.(i)) q) 5

let surface_merges () =
  exact "surface merge interleaved"
    (Sort.merge [| 1.0; 4.0; 7.0 |] [| 2.0; 3.0; 9.0 |])
    [| 1.0; 2.0; 3.0; 4.0; 7.0; 9.0 |] 6;
  exact "surface merge left empty" (Sort.merge [||] [| 1.0; 2.0 |])
    [| 1.0; 2.0 |] 2;
  exact "surface merge right empty" (Sort.merge [| 1.0; 2.0 |] [||])
    [| 1.0; 2.0 |] 2;
  let m = Sort.merge [| 1.0; 5.0 |] [| 2.0; 3.0; 4.0; 6.0 |] in
  check "surface merge length" (Array.length m = 6);
  ordered "surface merge order" true m 6

let surface_empty () =
  let a = [||] in
  Sort.sort a;
  Sort.sort ~desc:true ~carry:[||] a;
  Sort.permute a [||];
  check "surface empty index" (Array.length (Sort.sort_index [||]) = 0);
  check "surface empty merge" (Array.length (Sort.merge [||] [||]) = 0)

let surface_errors () =
  raises "surface sort short carry" (fun () ->
      Sort.sort ~carry:[| 1.0 |] [| 1.0; 2.0 |]);
  raises "surface permute wrong length" (fun () ->
      Sort.permute [| 1.0; 2.0 |] [| 0 |]);
  raises "surface permute not a permutation" (fun () ->
      Sort.permute [| 1.0; 2.0 |] [| 0; 0 |]);
  raises "surface permute out of range" (fun () ->
      Sort.permute [| 1.0; 2.0 |] [| 0; 2 |])

let () =
  one ();
  two_up ();
  two_down ();
  two_equal ();
  all_equal ();
  already_sorted ();
  reverse_sorted ();
  saw_tooth ();
  two_values ();
  organ_pipe ();
  signed_mix ();
  random_sizes ();
  tail_untouched ();
  dsort_ignores_dy ();
  dpperm_identity ();
  dpperm_cycle ();
  dpperm_swaps ();
  errors ();
  merge_plain ();
  merge_interleaved ();
  merge_disjoint_ranges ();
  merge_ties ();
  merge_singletons ();
  merge_empty ();
  merge_offsets ();
  merge_random ();
  surface_sorts ();
  surface_carries ();
  surface_indexes ();
  surface_merges ();
  surface_empty ();
  surface_errors ();
  if !bad = 0 then print_string "sort: all checks passed\n" else exit 1
