
let naif = 8
let naig = 8
let naip1 = 26
let naip2 = 15
let nbif = 8
let nbig = 8
let nbif2 = 10
let nbig2 = 9
let nbip1 = 24
let nbip2 = 32
let nam20 = 17
let nath0 = 15
let nam21 = 25
let nath1 = 23
let nam22 = 34
let nath2 = 33
let ntpsi = 23
let ntapsi = 6
let ncot = 8

let pi = 3.14159265358979323846264338327950
let pi4 = 0.78539816339744830961566084581988
let pi2rec = 0.011619772367581343075535053490057
let atr = 8.75069057084843450880771988210148
let btr = -2.09383632135605431360096498526268

let d9aimp_xsml = -1.0 /. (eps_2_dp ** 0.3333)

let dai_x3sml = eps_2_dp ** 0.3334

let dai_xmax =
  let xmaxt = (-1.5 *. log tiny_dp) ** 0.6667 in
  xmaxt -. (xmaxt *. log xmaxt /. ((4.0 *. sqrt xmaxt) +. 1.0)) -. 0.01

let airy_eta = 0.1 *. eps_2_dp
let airy_xbig = huge_dp ** 0.6666
let daie_x3sml = airy_eta ** single 0.3333
let daie_x32sml = 1.3104 *. (daie_x3sml *. daie_x3sml)
let bairy_x3sml = airy_eta ** 0.3333
let dbie_x32sml = 1.3104 *. (bairy_x3sml *. bairy_x3sml)
let dbi_xmax = (1.5 *. log huge_dp) ** 0.6666

let dcot_xmax = 1.0 /. eps_dp
let dcot_xsml = sqrt (3.0 *. eps_2_dp)
let dcot_xmin = exp (Float.max (log tiny_dp) (-.(log huge_dp)) +. 0.01)

let dpsi_xbig = 1.0 /. sqrt eps_2_dp

let gamln_con = 1.83787706640934548
let gamln_wdtol = Float.max eps_dp 0.5e-18

let gamln_zmin =
  let rln = log10_radix_dp *. 53.0 in
  let fln = Float.max (Float.min rln 20.0) 3.0 -. 3.0 in
  float_of_int (int_of_float (1.8 +. (0.3875 *. fln) +. 1.0))

let drf_errtol = (4.0 *. eps_2_dp) ** (1.0 /. 6.0)
let drf_lolim = 5.0 *. tiny_dp
let drf_uplim = huge_dp /. 5.0
let drf_c1 = 1.0 /. 24.0
let drf_c2 = 3.0 /. 44.0
let drf_c3 = 1.0 /. 14.0

let drd_errtol = (eps_2_dp /. 3.0) ** (1.0 /. 6.0)
let drd_lolim = 2.0 /. (huge_dp ** (2.0 /. 3.0))
let drd_uplim = (0.10 *. drd_errtol /. tiny_dp) ** (2.0 /. 3.0)
let drd_c1 = 3.0 /. 14.0
let drd_c2 = 1.0 /. 6.0
let drd_c3 = 9.0 /. 22.0
let drd_c4 = 3.0 /. 26.0

let drc_errtol = (eps_2_dp /. 16.0) ** (1.0 /. 6.0)
let drc_lolim = 5.0 *. tiny_dp
let drc_uplim = huge_dp /. 5.0
let drc_c1 = 1.0 /. 7.0
let drc_c2 = 9.0 /. 22.0

let carlson_maxit = 100

let d9aimp x =
  if x >= -1.0 then invalid_arg "d9aimp: x must be < -1";
  let ampl, theta =
    if x >= -2.0 then begin
      let z = ((16.0 /. cube x) +. 9.0) /. 7.0 in
      (0.3125 +. csevl z am22cs nam22, -0.625 +. csevl z ath2cs nath2)
    end else if x >= -4.0 then begin
      let z = ((128.0 /. cube x) +. 9.0) /. 7.0 in
      (0.3125 +. csevl z am21cs nam21, -0.625 +. csevl z ath1cs nath1)
    end else begin
      let z = if x > d9aimp_xsml then (128.0 /. cube x) +. 1.0 else 1.0 in
      (0.3125 +. csevl z am20cs nam20, -0.625 +. csevl z ath0cs nath0)
    end
  in
  let sqrtx = sqrt (-.x) in
  (sqrt (ampl /. sqrtx), pi4 -. (x *. sqrtx *. theta))

let daie x =
  if x < -1.0 then begin
    let xm, theta = d9aimp x in
    xm *. cos theta
  end else if x <= 1.0 then begin
    let z = if Float.abs x > daie_x3sml then cube x else 0.0 in
    let r =
      0.375 +. (csevl z aifcs naif -. (x *. (0.25 +. csevl z aigcs naig)))
    in
    if x > daie_x32sml then r *. exp (2.0 *. x *. sqrt x /. 3.0) else r
  end else if x <= 4.0 then begin
    let sqrtx = sqrt x in
    let z = ((16.0 /. (x *. sqrtx)) -. 9.0) /. 7.0 in
    (0.28125 +. csevl z aip1cs naip1) /. sqrt sqrtx
  end else begin
    let sqrtx = sqrt x in
    let z = if x < airy_xbig then (16.0 /. (x *. sqrtx)) -. 1.0 else -1.0 in
    (0.28125 +. csevl z aip2cs naip2) /. sqrt sqrtx
  end

let dai x =
  if x < -1.0 then begin
    let xm, theta = d9aimp x in
    xm *. cos theta
  end else if x <= 1.0 then begin
    let z = if Float.abs x > dai_x3sml then cube x else 0.0 in
    0.375 +. (csevl z aifcs naif -. (x *. (0.25 +. csevl z aigcs naig)))
  end else if x <= dai_xmax then daie x *. exp (-2.0 *. x *. sqrt x /. 3.0)
  else 0.0

let dbie x =
  if x < -1.0 then begin
    let xm, theta = d9aimp x in
    xm *. sin theta
  end else if x <= 1.0 then begin
    let z = if Float.abs x > bairy_x3sml then cube x else 0.0 in
    let r =
      0.625 +. csevl z bifcs nbif +. (x *. (0.4375 +. csevl z bigcs nbig))
    in
    if x > dbie_x32sml then r *. exp (-2.0 *. x *. sqrt x /. 3.0) else r
  end else if x <= 2.0 then begin
    let z = ((2.0 *. cube x) -. 9.0) /. 7.0 in
    exp (-2.0 *. x *. sqrt x /. 3.0)
    *. (1.125 +. csevl z bif2cs nbif2 +. (x *. (0.625 +. csevl z big2cs nbig2)))
  end else if x <= 4.0 then begin
    let sqrtx = sqrt x in
    let z = (atr /. (x *. sqrtx)) +. btr in
    (0.625 +. csevl z bip1cs nbip1) /. sqrt sqrtx
  end else begin
    let sqrtx = sqrt x in
    let z = if x < airy_xbig then (16.0 /. (x *. sqrtx)) -. 1.0 else -1.0 in
    (0.625 +. csevl z bip2cs nbip2) /. sqrt sqrtx
  end

let dbi x =
  if x < -1.0 then begin
    let xm, theta = d9aimp x in
    xm *. sin theta
  end else if x <= 1.0 then begin
    let z = if Float.abs x > bairy_x3sml then cube x else 0.0 in
    0.625 +. csevl z bifcs nbif +. (x *. (0.4375 +. csevl z bigcs nbig))
  end else if x <= 2.0 then begin
    let z = ((2.0 *. cube x) -. 9.0) /. 7.0 in
    1.125 +. csevl z bif2cs nbif2 +. (x *. (0.625 +. csevl z big2cs nbig2))
  end else if x <= dbi_xmax then dbie x *. exp (2.0 *. x *. sqrt x /. 3.0)
  else invalid_arg "dbi: x so big that bi overflows"

let dcot x =
  let ax = Float.abs x in
  if ax < dcot_xmin then
    invalid_arg "dcot: abs x is zero or so small dcot overflows";
  if ax > dcot_xmax then
    invalid_arg "dcot: no precision because abs x is too big";
  let ainty = Float.trunc ax in
  let yrem = ax -. ainty in
  let prodbg = 0.625 *. ainty in
  let ainty = Float.trunc prodbg in
  let y = (prodbg -. ainty) +. (0.625 *. yrem) +. (pi2rec *. ax) in
  let ainty2 = Float.trunc y in
  let ainty = ainty +. ainty2 in
  let y = y -. ainty2 in
  let ifn = int_of_float (Float.rem ainty 2.0) in
  let y = if ifn = 1 then 1.0 -. y else y in
  let r =
    if y <= dcot_xsml then 1.0 /. x
    else if y <= 0.25 then
      (0.5 +. csevl ((32.0 *. y *. y) -. 1.0) cotcs ncot) /. y
    else if y <= 0.5 then begin
      let c = (0.5 +. csevl ((8.0 *. y *. y) -. 1.0) cotcs ncot) /. (0.5 *. y) in
      ((c *. c) -. 1.0) *. 0.5 /. c
    end else begin
      let c = (0.5 +. csevl ((2.0 *. y *. y) -. 1.0) cotcs ncot) /. (0.25 *. y) in
      let c = ((c *. c) -. 1.0) *. 0.5 /. c in
      ((c *. c) -. 1.0) *. 0.5 /. c
    end
  in
  let r = Float.copy_sign r x in
  if ifn = 1 then -.r else r

let dpsi x =
  let ax = Float.abs x in
  if ax > 10.0 then begin
    let t = 10.0 /. ax in
    let aux =
      if ax < dpsi_xbig then csevl ((2.0 *. (t *. t)) -. 1.0) apsics ntapsi
      else 0.0
    in
    if x < 0.0 then
      log (Float.abs x) -. (0.5 /. x) +. aux -. (pi *. dcot (pi *. x))
    else log x -. (0.5 /. x) +. aux
  end else begin
    let n = int_of_float x in
    let n = if x < 0.0 then n - 1 else n in
    let y = x -. float_of_int n in
    let n = n - 1 in
    let s = csevl ((2.0 *. y) -. 1.0) psics ntpsi in
    if n = 0 then s
    else if n < 0 then begin
      let n = -n in
      if x = 0.0 then invalid_arg "dpsi: x is 0";
      if x < 0.0 && x +. float_of_int n -. 2.0 = 0.0 then
        invalid_arg "dpsi: x is a negative integer";
      let acc = ref s in
      for i = 1 to n do
        acc := !acc -. (1.0 /. (x +. float_of_int (i - 1)))
      done;
      !acc
    end else begin
      let acc = ref s in
      for i = 1 to n do
        acc := !acc +. (1.0 /. (y +. float_of_int i))
      done;
      !acc
    end
  end

let dgamln z =
  if z <= 0.0 then invalid_arg "dgamln: z <= 0";
  let nz = if z > 101.0 then 0 else int_of_float z in
  if z <= 101.0 && z -. float_of_int nz = 0.0 && nz <= 100 then gln.(nz - 1)
  else begin
    let zinc = if z >= gamln_zmin then 0.0 else gamln_zmin -. float_of_int nz in
    let zdmy = z +. zinc in
    let zp = 1.0 /. zdmy in
    let t1 = cf.(0) *. zp in
    let s = ref t1 in
    if zp >= gamln_wdtol then begin
      let zsq = zp *. zp in
      let tst = t1 *. gamln_wdtol in
      let zpk = ref zp and k = ref 1 and go = ref true in
      while !go && !k <= 21 do
        zpk := !zpk *. zsq;
        let trm = cf.(!k) *. !zpk in
        if Float.abs trm < tst then go := false
        else begin
          s := !s +. trm;
          incr k
        end
      done
    end;
    if zinc = 0.0 then begin
      let tlg = log z in
      (z *. (tlg -. 1.0)) +. (0.5 *. (gamln_con -. tlg)) +. !s
    end else begin
      let zp = ref 1.0 in
      for i = 1 to int_of_float zinc do
        zp := !zp *. (z +. float_of_int (i - 1))
      done;
      let tlg = log zdmy in
      (zdmy *. (tlg -. 1.0)) -. log !zp +. (0.5 *. (gamln_con -. tlg)) +. !s
    end
  end

let drf x y z =
  if Float.min x (Float.min y z) < 0.0 then invalid_arg "drf: min(x,y,z) < 0";
  if Float.max x (Float.max y z) > drf_uplim then
    invalid_arg "drf: max(x,y,z) > uplim";
  if Float.min (x +. y) (Float.min (x +. z) (y +. z)) < drf_lolim then
    invalid_arg "drf: min(x+y,x+z,y+z) < lolim";
  let xn = ref x and yn = ref y and zn = ref z in
  let res = ref 0.0 and converged = ref false and it = ref 0 in
  while (not !converged) && !it < carlson_maxit do
    incr it;
    let mu = (!xn +. !yn +. !zn) /. 3.0 in
    let xndev = 2.0 -. ((mu +. !xn) /. mu) in
    let yndev = 2.0 -. ((mu +. !yn) /. mu) in
    let zndev = 2.0 -. ((mu +. !zn) /. mu) in
    let epslon =
      Float.max (Float.abs xndev)
        (Float.max (Float.abs yndev) (Float.abs zndev))
    in
    if epslon < drf_errtol then begin
      let e2 = (xndev *. yndev) -. (zndev *. zndev) in
      let e3 = xndev *. yndev *. zndev in
      let s =
        1.0 +. (((drf_c1 *. e2) -. 0.10 -. (drf_c2 *. e3)) *. e2)
        +. (drf_c3 *. e3)
      in
      res := s /. sqrt mu;
      converged := true
    end else begin
      let xnroot = sqrt !xn and ynroot = sqrt !yn and znroot = sqrt !zn in
      let lamda = (xnroot *. (ynroot +. znroot)) +. (ynroot *. znroot) in
      xn := (!xn +. lamda) *. 0.250;
      yn := (!yn +. lamda) *. 0.250;
      zn := (!zn +. lamda) *. 0.250
    end
  done;
  if not !converged then invalid_arg "drf: iteration limit reached";
  !res

let drd x y z =
  if Float.min x y < 0.0 then invalid_arg "drd: min(x,y) < 0";
  if Float.max x (Float.max y z) > drd_uplim then
    invalid_arg "drd: max(x,y,z) > uplim";
  if Float.min (x +. y) z < drd_lolim then invalid_arg "drd: min(x+y,z) < lolim";
  let xn = ref x and yn = ref y and zn = ref z in
  let sigma = ref 0.0 and power4 = ref 1.0 in
  let res = ref 0.0 and converged = ref false and it = ref 0 in
  while (not !converged) && !it < carlson_maxit do
    incr it;
    let mu = (!xn +. !yn +. (3.0 *. !zn)) *. 0.20 in
    let xndev = (mu -. !xn) /. mu in
    let yndev = (mu -. !yn) /. mu in
    let zndev = (mu -. !zn) /. mu in
    let epslon =
      Float.max (Float.abs xndev)
        (Float.max (Float.abs yndev) (Float.abs zndev))
    in
    if epslon < drd_errtol then begin
      let ea = xndev *. yndev in
      let eb = zndev *. zndev in
      let ec = ea -. eb in
      let ed = ea -. (6.0 *. eb) in
      let ef = ed +. ec +. ec in
      let s1 =
        ed
        *. (-.drd_c1 +. (0.250 *. drd_c3 *. ed)
            -. (1.50 *. drd_c4 *. zndev *. ef))
      in
      let s2 =
        zndev
        *. ((drd_c2 *. ef)
            +. (zndev *. ((-.drd_c3 *. ec) +. (zndev *. drd_c4 *. ea))))
      in
      res :=
        (3.0 *. !sigma) +. (!power4 *. (1.0 +. s1 +. s2) /. (mu *. sqrt mu));
      converged := true
    end else begin
      let xnroot = sqrt !xn and ynroot = sqrt !yn and znroot = sqrt !zn in
      let lamda = (xnroot *. (ynroot +. znroot)) +. (ynroot *. znroot) in
      sigma := !sigma +. (!power4 /. (znroot *. (!zn +. lamda)));
      power4 := !power4 *. 0.250;
      xn := (!xn +. lamda) *. 0.250;
      yn := (!yn +. lamda) *. 0.250;
      zn := (!zn +. lamda) *. 0.250
    end
  done;
  if not !converged then invalid_arg "drd: iteration limit reached";
  !res

let drc x y =
  if x < 0.0 || y <= 0.0 then invalid_arg "drc: x < 0 or y <= 0";
  if Float.max x y > drc_uplim then invalid_arg "drc: max(x,y) > uplim";
  if x +. y < drc_lolim then invalid_arg "drc: x+y < lolim";
  let xn = ref x and yn = ref y in
  let res = ref 0.0 and converged = ref false and it = ref 0 in
  while (not !converged) && !it < carlson_maxit do
    incr it;
    let mu = (!xn +. !yn +. !yn) /. 3.0 in
    let sn = ((!yn +. mu) /. mu) -. 2.0 in
    if Float.abs sn < drc_errtol then begin
      let s =
        sn *. sn
        *. (0.30 +. (sn *. (drc_c1 +. (sn *. (0.3750 +. (sn *. drc_c2))))))
      in
      res := (1.0 +. s) /. sqrt mu;
      converged := true
    end else begin
      let lamda = (2.0 *. sqrt !xn *. sqrt !yn) +. !yn in
      xn := (!xn +. lamda) *. 0.250;
      yn := (!yn +. lamda) *. 0.250
    end
  done;
  if not !converged then invalid_arg "drc: iteration limit reached";
  !res
