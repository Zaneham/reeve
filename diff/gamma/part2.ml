let nti0 = 11
let ntai0 = 22
let ntai02 = 24
let nti1 = 11
let ntai1 = 22
let ntai12 = 24
let ntk0 = 10
let ntak0 = 17
let ntak02 = 14
let ntk1 = 10
let ntak1 = 18
let ntak12 = 14

let dcsevl x cs n =
  if Float.abs x > 1.0 then
    invalid_arg "dcsevl: x outside the interval (-1,+1)";
  let b0 = ref 0.0 and b1 = ref 0.0 and b2 = ref 0.0 in
  let twox = 2.0 *. x in
  for i = 1 to n do
    b2 := !b1;
    b1 := !b0;
    b0 := twox *. !b1 -. !b2 +. cs.(n - i)
  done;
  0.5 *. (!b0 -. !b2)

let gamln_con = 1.83787706640934548
let gamln_zmin =
  let rln = log10_radix_dp *. float_of_int digits_dp in
  let fln = Float.max (Float.min rln 20.0) 3.0 -. 3.0 in
  float_of_int (int_of_float (1.8 +. 0.3875 *. fln +. 1.0))

let dgamln z =
  if z <= 0.0 then invalid_arg "dgamln: z <= 0";
  let tabled =
    if z > 101.0 then None
    else
      let nz = int_of_float z in
      if z -. float_of_int nz > 0.0 || nz > 100 then None
      else Some gamln_tab.(nz - 1)
  in
  match tabled with
  | Some v -> v
  | None ->
    let wdtol = Float.max epsilon_float 0.5e-18 in
    let zinc = if z < gamln_zmin then gamln_zmin -. float_of_int (int_of_float z) else 0.0 in
    let zdmy = z +. zinc in
    let zp = ref (1.0 /. zdmy) in
    let t1 = cf_tab.(0) *. !zp in
    let s = ref t1 in
    if !zp >= wdtol then begin
      let zsq = !zp *. !zp in
      let tst = t1 *. wdtol in
      let k = ref 2 and brk = ref false in
      while not !brk && !k <= 22 do
        zp := !zp *. zsq;
        let trm = cf_tab.(!k - 1) *. !zp in
        if Float.abs trm < tst then brk := true
        else begin
          s := !s +. trm;
          incr k
        end
      done
    end;
    if zinc = 0.0 then
      let tlg = log z in
      (z *. (tlg -. 1.0)) +. (0.5 *. (gamln_con -. tlg)) +. !s
    else begin
      let prod = ref 1.0 in
      for i = 1 to int_of_float zinc do
        prod := !prod *. (z +. float_of_int (i - 1))
      done;
      let tlg = log zdmy in
      (zdmy *. (tlg -. 1.0)) -. log !prod +. (0.5 *. (gamln_con -. tlg)) +. !s
    end

let dgamma_pos a = exp (dgamln a)

let i0_xsml = sqrt (4.5 *. eps_2_dp)
let i0_xmax = log huge_dp
let i1_xsml = sqrt (4.5 *. eps_2_dp)
let i1_xmin = 2.0 *. tiny_dp
let i1_xmax = log huge_dp
let k_xsml = sqrt (4.0 *. eps_2_dp)
let k_xmaxt = -.log tiny_dp
let k_xmax = k_xmaxt -. (0.5 *. k_xmaxt *. log k_xmaxt /. (k_xmaxt +. 0.5))
let k1_xmin = exp (Float.max (log tiny_dp) (-.log huge_dp) +. 0.01)

let dbesi0e x =
  let y = Float.abs x in
  if y > 8.0 then (0.375 +. dcsevl (16.0 /. y -. 1.0) ai02cs ntai02) /. sqrt y
  else if y > 3.0 then
    (0.375 +. dcsevl ((48.0 /. y -. 11.0) /. 5.0) ai0cs ntai0) /. sqrt y
  else if y > i0_xsml then
    exp (-.y) *. (2.75 +. dcsevl (y *. y /. 4.5 -. 1.0) bi0cs nti0)
  else 1.0 -. x

let dbesi0 x =
  let y = Float.abs x in
  if y > i0_xmax then invalid_arg "dbesi0: abs x so big that i0 overflows"
  else if y > 3.0 then exp y *. dbesi0e x
  else if y > i0_xsml then 2.75 +. dcsevl (y *. y /. 4.5 -. 1.0) bi0cs nti0
  else 1.0

let dbesi1e x =
  let y = Float.abs x in
  if y > 3.0 then
    let v =
      if y <= 8.0 then (0.375 +. dcsevl ((48.0 /. y -. 11.0) /. 5.0) ai1cs ntai1) /. sqrt y
      else (0.375 +. dcsevl (16.0 /. y -. 1.0) ai12cs ntai12) /. sqrt y
    in
    Float.copy_sign v x
  else if x > i1_xmin then
    let v =
      if y > i1_xsml then x *. (0.875 +. dcsevl (y *. y /. 4.5 -. 1.0) bi1cs nti1)
      else 0.5 *. x
    in
    exp (-.y) *. v
  else 0.0

let dbesi1 x =
  let y = Float.abs x in
  if y > i1_xmax then invalid_arg "dbesi1: abs x so big that i1 overflows"
  else if y > 3.0 then exp y *. dbesi1e x
  else if y > i1_xsml then x *. (0.875 +. dcsevl (y *. y /. 4.5 -. 1.0) bi1cs nti1)
  else if y > i1_xmin then 0.5 *. x
  else 0.0

let dbesk0e x =
  if x <= 0.0 then invalid_arg "dbesk0e: x <= 0"
  else if x <= 2.0 then
    let y = if x > k_xsml then x *. x else 0.0 in
    exp x
    *. (-.log (0.5 *. x) *. dbesi0 x -. 0.25 +. dcsevl (0.5 *. y -. 1.0) bk0cs ntk0)
  else if x <= 8.0 then
    (1.25 +. dcsevl ((16.0 /. x -. 5.0) /. 3.0) ak0cs ntak0) /. sqrt x
  else (1.25 +. dcsevl (16.0 /. x -. 1.0) ak02cs ntak02) /. sqrt x

let dbesk0 x =
  if x <= 0.0 then invalid_arg "dbesk0: x is zero or negative"
  else if x <= 2.0 then
    let y = if x > k_xsml then x *. x else 0.0 in
    -.log (0.5 *. x) *. dbesi0 x -. 0.25 +. dcsevl (0.5 *. y -. 1.0) bk0cs ntk0
  else if x <= k_xmax then exp (-.x) *. dbesk0e x
  else 0.0

let dbesk1e x =
  if x <= 0.0 then invalid_arg "dbesk1e: x is zero or negative"
  else if x < k1_xmin then invalid_arg "dbesk1e: x so small that k1 overflows"
  else if x <= 2.0 then
    let y = if x > k_xsml then x *. x else 0.0 in
    exp x
    *. (log (0.5 *. x) *. dbesi1 x
        +. (0.75 +. dcsevl (0.5 *. y -. 1.0) bk1cs ntk1) /. x)
  else if x <= 8.0 then
    (1.25 +. dcsevl ((16.0 /. x -. 5.0) /. 3.0) ak1cs ntak1) /. sqrt x
  else (1.25 +. dcsevl (16.0 /. x -. 1.0) ak12cs ntak12) /. sqrt x

let dbesk1 x =
  if x <= 0.0 then invalid_arg "dbesk1: x is zero or negative"
  else if x < k1_xmin then invalid_arg "dbesk1: x so small that k1 overflows"
  else if x <= 2.0 then
    let y = if x > k_xsml then x *. x else 0.0 in
    log (0.5 *. x) *. dbesi1 x +. (0.75 +. dcsevl (0.5 *. y -. 1.0) bk1cs ntk1) /. x
  else if x <= k_xmax then exp (-.x) *. dbesk1e x
  else 0.0
