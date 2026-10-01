let jair_n1 = 14
let jair_n2 = 23
let jair_n3 = 19
let jair_n4 = 15
let jair_m1 = 12
let jair_m2 = 21
let jair_m3 = 17
let jair_m4 = 13
let jair_n1d = 14
let jair_n2d = 24
let jair_n3d = 19
let jair_n4d = 15
let jair_m1d = 12
let jair_m2d = 22
let jair_m3d = 17
let jair_m4d = 13

let jair_fpi12 = 1.30899693899575e+00
let jair_con2 = 5.03154716196777e+00
let jair_con3 = 3.80004589867293e-01
let jair_con4 = 8.33333333333333e-01
let jair_con5 = 8.66025403784439e-01

let clenshaw tab n m t tt =
  let j = ref n in
  let f1 = ref tab.(!j - 1) and f2 = ref 0.0 in
  for _i = 1 to m do
    decr j;
    let prev = !f1 in
    f1 := tt *. !f1 -. !f2 +. tab.(!j - 1);
    f2 := prev
  done;
  t *. !f1 -. !f2 +. tab.(0)

let djairy x rx c =
  let negative_tail () =
    let t = 0.4 *. c -. 1.0 in
    let tt = t +. t in
    let ai =
      clenshaw ajn_tab jair_n3 jair_m3 t tt -. x *. clenshaw ajp_tab jair_n3 jair_m3 t tt
    in
    let dai =
      x *. x *. clenshaw dajp_tab jair_n3d jair_m3d t tt
      +. clenshaw dajn_tab jair_n3d jair_m3d t tt
    in
    (ai, dai)
  in
  if x < 0.0 then
    if c > 5.0 then begin
      let t = 10.0 /. c -. 1.0 in
      let tt = t +. t in
      let temp1 = clenshaw a_tab jair_n4 jair_m4 t tt in
      let temp2 = clenshaw b_tab jair_n4 jair_m4 t tt in
      let rtrx = sqrt rx in
      let cv = c -. jair_fpi12 in
      let ccv = cos cv and scv = sin cv in
      let ai = (temp1 *. ccv -. temp2 *. scv) /. rtrx in
      let temp1 = clenshaw da_tab jair_n4d jair_m4d t tt in
      let temp2 = clenshaw db_tab jair_n4d jair_m4d t tt in
      let e1 = ccv *. jair_con5 +. 0.5 *. scv in
      let e2 = scv *. jair_con5 -. 0.5 *. ccv in
      (ai, (temp1 *. e1 -. temp2 *. e2) *. rtrx)
    end
    else negative_tail ()
  else if c > 5.0 then begin
    let t = 10.0 /. c -. 1.0 in
    let tt = t +. t in
    let rtrx = sqrt rx in
    let ec = exp (-.c) in
    let ai = ec *. clenshaw ak3_tab jair_n1 jair_m1 t tt /. rtrx in
    let dai = -.rtrx *. ec *. clenshaw dak3_tab jair_n1d jair_m1d t tt in
    (ai, dai)
  end
  else if x > 1.20 then begin
    let t = (x +. x -. jair_con2) *. jair_con3 in
    let tt = t +. t in
    let rtrx = sqrt rx in
    let ec = exp (-.c) in
    let ai = ec *. clenshaw ak2_tab jair_n2 jair_m2 t tt /. rtrx in
    let dai = -.ec *. clenshaw dak2_tab jair_n2d jair_m2d t tt *. rtrx in
    (ai, dai)
  end
  else begin
    let t = (x +. x -. 1.2) *. jair_con4 in
    let tt = t +. t in
    let ai = clenshaw ak1_tab jair_n1 jair_m1 t tt in
    let dai = -.clenshaw dak1_tab jair_n1d jair_m1d t tt in
    (ai, dai)
  end

let asyjy_tols = -6.90775527898214
let asyjy_con1 = 6.66666666666667e-01
let asyjy_con2 = 3.33333333333333e-01
let asyjy_con548 = 1.04166666666667e-01

let dasyjy funjy x fnu flgjy in_ y yoff wk =
  let tol = Float.max eps_2_dp 1.0e-15 in
  let elim =
    if flgjy = 1.0 then
      -2.303 *. (log10_radix_dp *. float_of_int min_exp_dp +. 3.0)
    else -2.303 *. log10_radix_dp *. float_of_int (min_exp_dp + digits_dp)
  in
  let fn = ref fnu in
  let iflw = ref 0 in
  let upol = Array.make 10 0.0 in
  let cr = Array.make 10 0.0 in
  let dr = Array.make 10 0.0 in
  let jn = ref 1 in
  while !iflw = 0 && !jn <= in_ do
    let xx = x /. !fn in
    wk.(0) <- 1.0 -. xx *. xx;
    let abw2 = Float.abs wk.(0) in
    wk.(1) <- sqrt abw2;
    wk.(6) <- !fn ** asyjy_con2;
    let outcome =
      if abw2 > 0.27750 then begin
        let tau = 1.0 /. wk.(1) in
        let t2 = 1.0 /. wk.(0) in
        if wk.(0) >= 0.0 then begin
          wk.(2) <- Float.abs (log ((1.0 +. wk.(1)) /. xx) -. wk.(1));
          wk.(3) <- wk.(2) *. !fn;
          let rcz = asyjy_con1 /. wk.(3) in
          if wk.(3) <= elim then begin
            let z32 = 1.5 *. wk.(2) in
            let rtz = z32 ** asyjy_con2 in
            wk.(6) <- !fn ** asyjy_con2;
            wk.(4) <- rtz *. wk.(6);
            wk.(5) <- wk.(4) *. wk.(4);
            `Block100 (tau, t2, rcz, rtz)
          end
          else `Overflow
        end
        else begin
          wk.(2) <- Float.abs (wk.(1) -. atan wk.(1));
          wk.(3) <- wk.(2) *. !fn;
          let rcz = -.asyjy_con1 /. wk.(3) in
          let z32 = 1.5 *. wk.(2) in
          let rtz = z32 ** asyjy_con2 in
          wk.(4) <- rtz *. wk.(6);
          wk.(5) <- -.wk.(4) *. wk.(4);
          `Block100 (tau, t2, rcz, rtz)
        end
      end
      else begin
        let sa = ref 0.0 in
        if abw2 <> 0.0 then sa := asyjy_tols /. log abw2;
        let sb = !sa in
        let kmax = Array.make 5 0 in
        for i = 0 to 4 do
          kmax.(i) <- int_of_float (Float.max !sa 2.0);
          sa := !sa +. sb
        done;
        let kb = ref kmax.(4) in
        let klast = !kb - 1 in
        let sg = ref gama_tab.(!kb - 1) in
        for _k = 1 to klast do
          decr kb;
          sg := !sg *. wk.(0) +. gama_tab.(!kb - 1)
        done;
        let z = wk.(0) *. !sg in
        let az = Float.abs z in
        let rtz = sqrt az in
        wk.(2) <- asyjy_con1 *. az *. rtz;
        wk.(3) <- wk.(2) *. !fn;
        wk.(4) <- rtz *. wk.(6);
        wk.(5) <- -.wk.(4) *. wk.(4);
        if z > 0.0 && wk.(3) > elim then `Overflow
        else begin
          if z > 0.0 then wk.(5) <- -.wk.(5);
          let phi = sqrt (sqrt (!sg +. !sg +. !sg +. !sg)) in
          let kb = ref kmax.(4) in
          let klast = !kb - 1 in
          let sb0 = ref beta_tab.(0).(!kb - 1) in
          for _k = 1 to klast do
            decr kb;
            sb0 := !sb0 *. wk.(0) +. beta_tab.(0).(!kb - 1)
          done;
          let rfn2 = 1.0 /. (!fn *. !fn) in
          let rden = ref 1.0 in
          let asum = ref 1.0 in
          let relb = tol *. Float.abs !sb0 in
          let bsum = ref !sb0 in
          let ksp1 = ref 1 in
          let ks = ref 1 and brk = ref false in
          while not !brk && !ks <= 4 do
            incr ksp1;
            rden := !rden *. rfn2;
            let kstemp = 5 - !ks in
            let kb = ref kmax.(kstemp - 1) in
            let klast = !kb - 1 in
            let sa = ref alfa_tab.(!ks - 1).(!kb - 1) in
            let sb = ref beta_tab.(!ksp1 - 1).(!kb - 1) in
            for _k = 1 to klast do
              decr kb;
              sa := !sa *. wk.(0) +. alfa_tab.(!ks - 1).(!kb - 1);
              sb := !sb *. wk.(0) +. beta_tab.(!ksp1 - 1).(!kb - 1)
            done;
            let ta = !sa *. !rden and tb = !sb *. !rden in
            asum := !asum +. ta;
            bsum := !bsum +. tb;
            if Float.abs ta <= tol && Float.abs tb <= relb then brk := true;
            incr ks
          done;
          `Ready (phi, !asum, !bsum /. (!fn *. wk.(6)))
        end
      end
    in
    let outcome =
      match outcome with
      | `Overflow -> `Overflow
      | `Ready (phi, asum, bsum) -> `Ready (phi, asum, bsum)
      | `Block100 (tau, t2, rcz, rtz) ->
        let phi = sqrt ((rtz +. rtz) *. tau) in
        upol.(0) <- 1.0;
        let tb = ref 1.0 in
        let asum = ref 1.0 in
        let tfn = tau /. !fn in
        let rfn2 = (1.0 /. !fn) *. (1.0 /. !fn) in
        let rden = ref 1.0 in
        upol.(1) <- (c_tab.(0) *. t2 +. c_tab.(1)) *. tfn;
        let crz32 = asyjy_con548 *. rcz in
        let bsum = ref (upol.(1) +. crz32) in
        let relb = tol *. Float.abs !bsum in
        let ap = ref tfn in
        let ks = ref 0 in
        let kp1 = ref 2 in
        let rzden = ref rcz in
        let l = ref 2 in
        let iseta = ref false and isetb = ref false in
        let lr = ref 2 and brk = ref false in
        while not !brk && !lr <= 8 do
          let lrp1 = !lr + 1 in
          for _k = !lr to lrp1 do
            incr ks;
            incr kp1;
            incr l;
            let s1 = ref c_tab.(!l - 1) in
            for _j = 2 to !kp1 do
              incr l;
              s1 := !s1 *. t2 +. c_tab.(!l - 1)
            done;
            ap := !ap *. tfn;
            upol.(!kp1 - 1) <- !ap *. !s1;
            cr.(!ks - 1) <- br_tab.(!ks - 1) *. !rzden;
            rzden := !rzden *. rcz;
            dr.(!ks - 1) <- ar_tab.(!ks - 1) *. !rzden
          done;
          let suma = ref upol.(lrp1 - 1) in
          let sumb = ref (upol.(!lr + 1) +. upol.(lrp1 - 1) *. crz32) in
          let ju = ref lrp1 in
          for jr = 1 to !lr do
            decr ju;
            suma := !suma +. cr.(jr - 1) *. upol.(!ju - 1);
            sumb := !sumb +. dr.(jr - 1) *. upol.(!ju - 1)
          done;
          rden := !rden *. rfn2;
          tb := -. !tb;
          if wk.(0) > 0.0 then tb := Float.abs !tb;
          if !rden < tol then begin
            if not !iseta then begin
              if Float.abs !suma < tol then iseta := true;
              asum := !asum +. !suma *. !tb
            end;
            if not !isetb then begin
              if Float.abs !sumb < relb then isetb := true;
              bsum := !bsum +. !sumb *. !tb
            end;
            if !iseta && !isetb then brk := true
          end
          else begin
            asum := !asum +. !suma *. !tb;
            bsum := !bsum +. !sumb *. !tb
          end;
          lr := !lr + 2
        done;
        let tbf = if wk.(0) > 0.0 then -.wk.(4) else wk.(4) in
        `Ready (phi, !asum, !bsum /. tbf)
    in
    match outcome with
    | `Overflow -> iflw := 1
    | `Ready (phi, asum, bsum) ->
      let fi, dfi = funjy wk.(5) wk.(4) wk.(3) in
      let ta = 1.0 /. tol in
      let tb = tiny_dp *. ta *. 1.0e+3 in
      let fi, dfi, phi =
        if Float.abs fi <= tb then (fi *. ta, dfi *. ta, phi *. tol) else (fi, dfi, phi)
      in
      y.(yoff + !jn - 1) <- flgjy *. phi *. (fi *. asum +. dfi *. bsum) /. wk.(6);
      fn := !fn -. flgjy;
      incr jn
  done;
  !iflw

let dasyik x fnu kode flgik ra arg in_ y yoff =
  let tol = Float.max eps_2_dp 1.0e-15 in
  let fn = ref fnu in
  let kk = int_of_float ((3.0 -. flgik) /. 2.0) in
  for jn = 1 to in_ do
    if jn <> 1 then begin
      fn := !fn -. flgik;
      let z = x /. !fn in
      ra := sqrt (1.0 +. z *. z);
      let gln = log ((1.0 +. !ra) /. z) in
      let etx = float_of_int (kode - 1) in
      let t = !ra *. (1.0 -. etx) +. etx /. (z +. !ra) in
      arg := !fn *. (t -. gln) *. flgik
    end;
    let coef = exp !arg in
    let t0 = 1.0 /. !ra in
    let t2 = t0 *. t0 in
    let t = Float.copy_sign (t0 /. !fn) flgik in
    let s2 = ref 1.0 in
    let ap = ref 1.0 in
    let l = ref 0 in
    let k = ref 2 and brk = ref false in
    while not !brk && !k <= 11 do
      incr l;
      let s1 = ref c_tab.(!l - 1) in
      for _j = 2 to !k do
        incr l;
        s1 := !s1 *. t2 +. c_tab.(!l - 1)
      done;
      ap := !ap *. t;
      let ak = !ap *. !s1 in
      s2 := !s2 +. ak;
      if Float.max (Float.abs ak) (Float.abs !ap) < tol then brk := true;
      incr k
    done;
    y.(yoff + jn - 1) <- !s2 *. coef *. sqrt (Float.abs t) *. con_ik.(kk - 1)
  done
