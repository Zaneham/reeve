(* Reeve, Copyright 2026 Zane Hambly.
   Extended-range arithmetic of Smith, Olver and Lozier, and the Legendre
   functions of the first and second kind built on top of it, ported from the
   SLATEC double precision routines. *)

open Mach

type xnum = { x : float; ix : int }

type env = {
  nbitsf : int;
  l : int;
  l2 : int;
  kmax : int;
  nlg102 : int;
  mlg102 : int;
  lg102 : int array;
  log2r : int;
  radix : float;
  radixl : float;
  rad2l : float;
  dlg10r : float;
}

type kind = Pneg | Q | Ppos | Pnorm
type mode = By_x | By_theta

let zero = { x = 0.0; ix = 0 }

(* ---- Extended-range arithmetic ---- *)

let log102 =
  [|
    301; 029; 995; 663; 981; 195; 213; 738; 894; 724; 493; 026; 768; 189; 881;
    462; 108; 541; 310; 428;
  |]

let dxset irad nradpl dzero nbits =
  let iradx = if irad = 0 then 2 else irad in
  let nrdplc = if nradpl = 0 then 53 else nradpl in
  let iminex = if dzero = 0.0 then -1021 else 0 in
  let imaxex = if dzero = 0.0 then 1024 else 0 in
  let nbitsx = if nbits = 0 then Sys.int_size - 1 else nbits in
  let log2r =
    if iradx = 2 then 1
    else if iradx = 4 then 2
    else if iradx = 8 then 3
    else if iradx = 16 then 4
    else invalid_arg "dxset: improper value of irad, argument 1"
  in
  let nbitsf = log2r * nrdplc in
  let radix = float_of_int iradx in
  let dlg10r = log10 radix in
  let lx =
    if dzero <> 0.0 then int_of_float (0.5 *. log10 dzero /. dlg10r) - 1
    else min ((1 - iminex) / 2) ((imaxex - 1) / 2)
  in
  let l2 = 2 * lx in
  if lx < 4 then invalid_arg "dxset: improper value of dzero, argument 3";
  let l = lx in
  let radixl = ldexp 1.0 (log2r * l) in
  let rad2l = radixl *. radixl in
  if rad2l = infinity then
    invalid_arg "dxset: radix**(2*l) overflows, arguments 1 and 3";
  if nbitsx < 15 || nbitsx > Sys.int_size - 1 then
    invalid_arg "dxset: improper value of nbits, argument 4";
  let kmax = (1 lsl (nbitsx - 1)) - l2 in
  let nb = (nbitsx - 1) / 2 in
  let mlg102 = 1 lsl nb in
  if nrdplc * log2r < 1 || nrdplc * log2r > 120 then
    invalid_arg "dxset: improper value of nradpl, argument 2";
  let nlg102 = (nrdplc * log2r / nb) + 3 in
  let np1 = nlg102 + 1 in
  let lgtemp = Array.make 20 0 in
  let ic = ref 0 in
  for ii = 1 to 20 do
    let i = 21 - ii in
    let it = (log2r * log102.(i - 1)) + !ic in
    ic := it / 1000;
    lgtemp.(i - 1) <- it mod 1000
  done;
  let lg102 = Array.make 21 0 in
  lg102.(0) <- !ic;
  for i = 2 to np1 do
    let lg102x = ref 0 in
    for _j = 1 to nb do
      ic := 0;
      for kk = 1 to 20 do
        let k = 21 - kk in
        let it = (2 * lgtemp.(k - 1)) + !ic in
        ic := it / 1000;
        lgtemp.(k - 1) <- it mod 1000
      done;
      lg102x := (2 * !lg102x) + !ic
    done;
    lg102.(i - 1) <- !lg102x
  done;
  if nrdplc >= l then invalid_arg "dxset: nradpl >= l, argument 2";
  if 6 * l > kmax then invalid_arg "dxset: 6*l > kmax, argument 4";
  {
    nbitsf;
    l;
    l2;
    kmax;
    nlg102;
    mlg102;
    lg102;
    log2r;
    radix;
    radixl;
    rad2l;
    dlg10r;
  }

let radpow env k = ldexp 1.0 (env.log2r * k)

(* Zero keeps its sign, since DXADJ and DXRED both clear IX without touching
   X. *)
let dxadj env a =
  let overflow () = invalid_arg "dxadj: overflow in auxiliary index" in
  if a.x = 0.0 then { a with ix = 0 }
  else if Float.abs a.x >= 1.0 then
    if Float.abs a.x < env.radixl then
      if abs a.ix > env.kmax then overflow () else a
    else if a.ix > 0 && a.ix > env.kmax - env.l2 then overflow ()
    else { x = a.x /. env.rad2l; ix = a.ix + env.l2 }
  else if env.radixl *. Float.abs a.x >= 1.0 then
    if abs a.ix > env.kmax then overflow () else a
  else if a.ix < 0 && a.ix < -env.kmax + env.l2 then overflow ()
  else { x = a.x *. env.rad2l; ix = a.ix - env.l2 }

let dxadd env a b =
  let l = env.l and radixl = env.radixl and rad2l = env.rad2l in
  let x = a.x and ix = a.ix and y = b.x and iy = b.ix in
  let res =
    if x = 0.0 then { x = y; ix = iy }
    else if y = 0.0 then a
    else if
      (not ((ix >= 0 && iy >= 0) || (ix < 0 && iy < 0)))
      && (abs ix > 6 * l || abs iy > 6 * l)
    then if ix >= 0 then a else { x = y; ix = iy }
    else
      let i = ix - iy in
      if i = 0 then
        if Float.abs x > 1.0 && Float.abs y > 1.0 then
          { x = (x /. radixl) +. (y /. radixl); ix = ix + l }
        else if Float.abs x < 1.0 && Float.abs y < 1.0 then
          { x = (x *. radixl) +. (y *. radixl); ix = ix - l }
        else { x = x +. y; ix }
      else begin
        let s, is, t = if i < 0 then (y, iy, x) else (x, ix, y) in
        let i1 = abs i / l and i2 = abs i mod l in
        let j, t =
          if Float.abs t >= radixl then
            if i1 - 2 >= 0 then (i1 - 2, t *. radpow env (-i2) /. rad2l)
            else if i1 - 1 >= 0 then (i1 - 1, t *. radpow env (-i2) /. radixl)
            else (i1, t *. radpow env (-i2))
          else if Float.abs t >= 1.0 then
            if i1 - 1 >= 0 then (i1 - 1, t *. radpow env (-i2) /. radixl)
            else (i1, t *. radpow env (-i2))
          else if radixl *. Float.abs t >= 1.0 then (i1, t *. radpow env (-i2))
          else (i1 + 1, t *. radpow env (l - i2))
        in
        let shifted =
          if j = 0 then Some s
          else if Float.abs s >= radixl || j > 3 then None
          else begin
            let nmul =
              if Float.abs s >= 1.0 then if j = 1 then 1 else 0
              else if radixl *. Float.abs s >= 1.0 then
                if j = 1 then 1 else if j = 2 then 2 else 0
              else j
            in
            if nmul = 0 then None
            else begin
              let sv = ref s in
              for _n = 1 to nmul do
                sv := !sv *. radixl
              done;
              Some !sv
            end
          end
        in
        match shifted with
        | None -> { x = s; ix = is }
        | Some s ->
          if Float.abs s > 1.0 && Float.abs t > 1.0 then
            { x = (s /. radixl) +. (t /. radixl); ix = is - (j * l) + l }
          else if Float.abs s < 1.0 && Float.abs t < 1.0 then
            { x = (s *. radixl) +. (t *. radixl); ix = is - (j * l) - l }
          else { x = s +. t; ix = is - (j * l) }
      end
  in
  dxadj env res

let dxred env a =
  if a.x = 0.0 then { a with ix = 0 }
  else begin
    let xa = ref (Float.abs a.x) in
    let keep = ref true in
    if a.ix <> 0 then begin
      let ixa = abs a.ix in
      let ixa1 = ref (ixa / env.l2) and ixa2 = ixa mod env.l2 in
      let bound = (2098 / env.l2) + 3 in
      let n = ref 0 in
      if a.ix > 0 then begin
        while !xa >= 1.0 && !n < bound do
          xa := !xa /. env.rad2l;
          incr ixa1;
          incr n
        done;
        if !xa >= 1.0 then invalid_arg "dxred: scaling loop hit its bound";
        xa := !xa *. radpow env ixa2;
        let i = ref 1 in
        while !keep && !i <= !ixa1 do
          if !xa > 1.0 then keep := false
          else begin
            xa := !xa *. env.rad2l;
            incr i
          end
        done
      end
      else begin
        while !xa <= 1.0 && !n < bound do
          xa := !xa *. env.rad2l;
          incr ixa1;
          incr n
        done;
        if !xa <= 1.0 then invalid_arg "dxred: scaling loop hit its bound";
        xa := !xa /. radpow env ixa2;
        let i = ref 1 in
        while !keep && !i <= !ixa1 do
          if !xa < 1.0 then keep := false
          else begin
            xa := !xa /. env.rad2l;
            incr i
          end
        done
      end
    end;
    if !keep && !xa <= env.rad2l && (!xa > 1.0 || env.rad2l *. !xa >= 1.0) then
      { x = Float.copy_sign !xa a.x; ix = 0 }
    else a
  end

let dxc210 env k =
  if k = 0 then (1.0, 0)
  else begin
    let m = env.mlg102 in
    let ka = abs k in
    let ka1 = ka / m and ka2 = ka mod m in
    if ka1 >= m then invalid_arg "dxc210: k too large";
    let nm1 = env.nlg102 - 1 and np1 = env.nlg102 + 1 in
    let it = ka2 * env.lg102.(np1 - 1) in
    let ic = ref (it / m) in
    let z = ref (float_of_int (it mod m)) in
    let fm = float_of_int m in
    let ja = ref 0 in
    if ka1 > 0 then begin
      for ii = 1 to nm1 do
        let i = np1 - ii in
        let it = (ka2 * env.lg102.(i - 1)) + (ka1 * env.lg102.(i)) + !ic in
        ic := it / m;
        z := (!z /. fm) +. float_of_int (it mod m)
      done;
      ja := (ka * env.lg102.(0)) + (ka1 * env.lg102.(1)) + !ic
    end
    else begin
      for ii = 1 to nm1 do
        let i = np1 - ii in
        let it = (ka2 * env.lg102.(i - 1)) + !ic in
        ic := it / m;
        z := (!z /. fm) +. float_of_int (it mod m)
      done;
      ja := (ka * env.lg102.(0)) + !ic
    end;
    let z = !z /. fm in
    if k > 0 then (10.0 ** (z -. 1.0), !ja + 1) else (10.0 ** -.z, - !ja)
  end

let ispace = 1

let dxcon env a =
  let a = dxred env a in
  if a.ix = 0 then a
  else begin
    let a = dxadj env a in
    let case2 = a.ix >= 0 in
    let xv = ref a.x and ixv = ref a.ix in
    if case2 then begin
      if Float.abs !xv < 1.0 then begin
        xv := !xv *. env.radixl;
        ixv := !ixv - env.l
      end
    end
    else if Float.abs !xv >= 1.0 then begin
      xv := !xv /. env.radixl;
      ixv := !ixv + env.l
    end;
    let maxit = env.l + 64 in
    let i = ref (int_of_float (log10 (Float.abs !xv) /. env.dlg10r)) in
    let av = ref (radpow env !i) in
    let n = ref 0 in
    if case2 then begin
      while !av > Float.abs !xv && !n < maxit do
        decr i;
        av := !av /. env.radix;
        incr n
      done;
      while Float.abs !xv >= env.radix *. !av && !n < maxit do
        incr i;
        av := !av *. env.radix;
        incr n
      done
    end
    else begin
      while !av > env.radix *. Float.abs !xv && !n < maxit do
        decr i;
        av := !av /. env.radix;
        incr n
      done;
      while Float.abs !xv >= !av && !n < maxit do
        incr i;
        av := !av *. env.radix;
        incr n
      done
    end;
    if !n >= maxit then invalid_arg "dxcon: radix scaling hit its bound";
    let itemp = ref (int_of_float (float_of_int ispace /. env.dlg10r)) in
    let aw = ref (radpow env !itemp) in
    let b = 10.0 ** float_of_int ispace in
    let n = ref 0 in
    while !aw > b && !n < maxit do
      decr itemp;
      aw := !aw /. env.radix;
      incr n
    done;
    while b >= !aw *. env.radix && !n < maxit do
      incr itemp;
      aw := !aw *. env.radix;
      incr n
    done;
    if !n >= maxit then invalid_arg "dxcon: decimal scaling hit its bound";
    if !itemp > 0 then begin
      let i1 = !i / !itemp in
      xv := !xv *. radpow env (-(i1 * !itemp));
      ixv := !ixv + (i1 * !itemp);
      let zz, j = dxc210 env !ixv in
      let j1 = j / ispace in
      let j2 = j - (j1 * ispace) in
      xv := !xv *. zz *. (10.0 ** float_of_int j2);
      ixv := j1 * ispace
    end
    else begin
      xv := !xv *. radpow env (- !i);
      ixv := !ixv + !i;
      let zz, j = dxc210 env !ixv in
      xv := !xv *. zz;
      ixv := j
    end;
    let n = ref 0 in
    if case2 then
      while 10.0 *. Float.abs !xv >= b && !n < maxit do
        xv := !xv /. b;
        ixv := !ixv + ispace;
        incr n
      done
    else
      while b *. Float.abs !xv < 1.0 && !n < maxit do
        xv := !xv *. b;
        ixv := !ixv - ispace;
        incr n
      done;
    if !n >= maxit then
      invalid_arg "dxcon: decimal normalisation hit its bound";
    { x = !xv; ix = !ixv }
  end

(* ---- Psi for the Legendre series ---- *)

let xpsi_cnum =
  [|
    1.0; -1.0; 1.0; -1.0; 1.0; -691.0; 1.0; -3617.0; 43867.0; -174611.0;
    77683.0; -236364091.0;
  |]

let xpsi_cdenom =
  [|
    12.0; 120.0; 252.0; 240.0; 132.0; 32760.0; 12.0; 8160.0; 14364.0; 6600.0;
    276.0; 65520.0;
  |]

let dxpsi a ipsik ipsix =
  let n = max 0 (ipsix - int_of_float a) in
  let b = float_of_int n +. a in
  let c = ref 0.0 in
  for i = 1 to ipsik - 1 do
    let k = ipsik - i in
    c := (!c +. (xpsi_cnum.(k - 1) /. xpsi_cdenom.(k - 1))) /. (b *. b)
  done;
  let v = log b -. (!c +. (0.5 /. b)) in
  if n = 0 then v
  else begin
    let s = ref 0.0 in
    for m = 1 to n do
      s := !s +. (1.0 /. (float_of_int (n - m) +. a))
    done;
    v -. !s
  end

(* ---- Legendre functions ---- *)

let dxpqnu env nu1 nu2 mu0 theta kind pqa =
  let j0 = env.nbitsf in
  let ipsik = 1 + (env.nbitsf / 10) in
  let ipsix = 5 * ipsik in
  let nu = ref (Float.rem nu1 1.0) in
  if !nu >= 0.5 then nu := !nu -. 1.0;
  (match kind with
  | Q -> ()
  | Pneg | Ppos | Pnorm -> if !nu > -0.5 then nu := !nu -. 1.0);
  let mu = ref mu0 in
  let dmu = ref (float_of_int mu0) in
  let factmu = ref 1.0 and ifa = ref 0 in
  for i = 1 to mu0 do
    factmu := !factmu *. float_of_int i;
    let r = dxadj env { x = !factmu; ix = !ifa } in
    factmu := r.x;
    ifa := r.ix
  done;
  let x = cos theta in
  let sh = sin (theta /. 2.0) in
  let y = sh *. sh in
  let rt = tan (theta /. 2.0) in
  let pq = ref zero and pq2 = ref zero in
  for j = 1 to 2 do
    (match kind with
    | Pneg | Ppos | Pnorm ->
      let pqx = ref 1.0 and ipq = ref 0 in
      let a = ref 1.0 and ia = ref 0 in
      let i = ref 2 and stop = ref false in
      while (not !stop) && !i <= j0 do
        let di = float_of_int !i in
        a :=
          !a *. y *. (di -. 2.0 -. !nu) *. (di -. 1.0 +. !nu)
          /. ((di -. 1.0 +. !dmu) *. (di -. 1.0));
        let r = dxadj env { x = !a; ix = !ia } in
        a := r.x;
        ia := r.ix;
        if !a = 0.0 then stop := true
        else begin
          let s = dxadd env { x = !pqx; ix = !ipq } { x = !a; ix = !ia } in
          pqx := s.x;
          ipq := s.ix;
          incr i
        end
      done;
      if !mu > 0 then begin
        let x1 = ref !pqx in
        for _i = 1 to !mu do
          x1 := !x1 *. rt;
          let r = dxadj env { x = !x1; ix = !ipq } in
          x1 := r.x;
          ipq := r.ix
        done;
        let r = dxadj env { x = !x1 /. !factmu; ix = !ipq - !ifa } in
        pqx := r.x;
        ipq := r.ix
      end;
      pq := { x = !pqx; ix = !ipq }
    | Q ->
      let z = -.log rt in
      let w = dxpsi (!nu +. 1.0) ipsik ipsix in
      let xs = 1.0 /. sin theta in
      let pqx = ref 0.0 and ipq = ref 0 in
      let a = ref 1.0 and ia = ref 0 in
      for k = 1 to j0 do
        let flok = float_of_int k in
        if k <> 1 then begin
          a :=
            !a *. y *. (flok -. 2.0 -. !nu) *. (flok -. 1.0 +. !nu)
            /. ((flok -. 1.0 +. !dmu) *. (flok -. 1.0));
          let r = dxadj env { x = !a; ix = !ia } in
          a := r.x;
          ia := r.ix
        end;
        let x1 =
          if !mu >= 1 then
            ((!nu *. (!nu +. 1.0) *. (z -. w +. dxpsi flok ipsik ipsix))
            +. ((!nu -. flok +. 1.0) *. (!nu +. flok) /. (2.0 *. flok)))
            *. !a
          else (dxpsi flok ipsik ipsix -. w +. z) *. !a
        in
        let s = dxadd env { x = !pqx; ix = !ipq } { x = x1; ix = !ia } in
        pqx := s.x;
        ipq := s.ix
      done;
      if !mu >= 1 then begin
        pqx := -.rt *. !pqx;
        let s = dxadd env { x = !pqx; ix = !ipq } { x = -.xs; ix = 0 } in
        pqx := s.x;
        ipq := s.ix
      end;
      if j = 2 then begin
        mu := - !mu;
        dmu := -. !dmu
      end;
      pq := { x = !pqx; ix = !ipq });
    if j = 1 then pq2 := !pq;
    nu := !nu +. 1.0
  done;
  let kk = ref 0 and fin = ref false in
  if not (!nu -. 1.5 < nu1) then begin
    incr kk;
    pqa.(!kk - 1) <- !pq2;
    if !nu > nu2 +. 0.5 then fin := true
  end;
  let pq1 = ref zero in
  let bound = int_of_float (Float.max 0.0 (nu2 -. !nu)) + 8 in
  let it = ref 0 in
  while (not !fin) && !it <= bound do
    incr it;
    pq1 := !pq;
    if not (!nu < nu1 +. 0.5) then begin
      incr kk;
      pqa.(!kk - 1) <- !pq;
      if !nu > nu2 +. 0.5 then fin := true
    end;
    if not !fin then begin
      let x1 = (2.0 *. !nu -. 1.0) /. (!nu +. !dmu) *. x *. !pq1.x in
      let x2 = (!nu -. 1.0 -. !dmu) /. (!nu +. !dmu) *. !pq2.x in
      let s = dxadd env { x = x1; ix = !pq1.ix } { x = -.x2; ix = !pq2.ix } in
      pq := dxadj env s;
      nu := !nu +. 1.0;
      pq2 := !pq1
    end
  done;
  if not !fin then invalid_arg "dxpqnu: nu recurrence hit its bound"

let dxpmu env nu1 nu2 mu1 mu2 theta x sx kind pqa =
  dxpqnu env nu1 nu2 mu2 theta kind pqa;
  let p0 = pqa.(0) in
  dxpqnu env nu1 nu2 (mu2 - 1) theta kind pqa;
  let n = mu2 - mu1 + 1 in
  let p1 = pqa.(0) in
  pqa.(n - 1) <- p0;
  if n > 1 then begin
    pqa.(n - 2) <- p1;
    for jf = n - 2 downto 1 do
      let dmu = float_of_int (mu1 + jf) in
      let x1 = 2.0 *. dmu *. x *. sx *. pqa.(jf).x in
      let x2 = -.(nu1 -. dmu) *. (nu1 +. dmu +. 1.0) *. pqa.(jf + 1).x in
      let s =
        dxadd env { x = x1; ix = pqa.(jf).ix } { x = x2; ix = pqa.(jf + 1).ix }
      in
      pqa.(jf - 1) <- dxadj env s
    done
  end

let dxpmup env nu1 nu2 mu1 mu2 pqa =
  let nu = ref nu1 in
  let mu = ref mu1 in
  let dmu = ref (float_of_int mu1) in
  let n = int_of_float (nu2 -. nu1 +. 0.1) + (mu2 - mu1) + 1 in
  let j = ref 1 in
  let early = ref false in
  if Float.rem (single nu1) 1.0 = 0.0 then begin
    let going = ref true in
    while !going && not !early do
      if !dmu < !nu +. 1.0 then going := false
      else begin
        pqa.(!j - 1) <- zero;
        incr j;
        if !j > n then early := true
        else begin
          if nu2 -. nu1 > 0.5 then nu := !nu +. 1.0;
          if mu2 > mu1 then mu := !mu + 1
        end
      end
    done
  end;
  if not !early then begin
    let prod = ref 1.0 and iprod = ref 0 in
    for lf = 1 to 2 * !mu do
      prod := !prod *. (!dmu -. !nu -. float_of_int lf);
      let r = dxadj env { x = !prod; ix = !iprod } in
      prod := r.x;
      iprod := r.ix
    done;
    for i = !j to n do
      if !mu <> 0 then begin
        let sgn = if !mu land 1 = 0 then 1.0 else -1.0 in
        let v = pqa.(i - 1) in
        pqa.(i - 1) <- dxadj env { x = v.x *. !prod *. sgn; ix = v.ix + !iprod }
      end;
      if nu2 -. nu1 > 0.5 then begin
        prod := !prod *. (-. !dmu -. !nu -. 1.0) /. (!dmu -. !nu -. 1.0);
        let r = dxadj env { x = !prod; ix = !iprod } in
        prod := r.x;
        iprod := r.ix;
        nu := !nu +. 1.0
      end
      else begin
        prod := (!dmu -. !nu) *. !prod *. (-. !dmu -. !nu -. 1.0);
        let r = dxadj env { x = !prod; ix = !iprod } in
        prod := r.x;
        iprod := r.ix;
        mu := !mu + 1;
        dmu := !dmu +. 1.0
      end
    done
  end

let dxpnrm env nu1 nu2 mu1 mu2 pqa =
  let l = int_of_float (float_of_int (mu2 - mu1) +. (nu2 -. nu1 +. 1.5)) in
  let mu = ref mu1 in
  let dmu = ref (float_of_int mu1) in
  let nu = ref nu1 in
  let j = ref 1 in
  let early = ref false in
  let going = ref true in
  while !going && not !early do
    if !dmu <= !nu then going := false
    else begin
      pqa.(!j - 1) <- zero;
      incr j;
      if !j > l then early := true
      else begin
        if mu2 > mu1 then dmu := !dmu +. 1.0;
        if nu2 -. nu1 > 0.5 then nu := !nu +. 1.0
      end
    end
  done;
  if not !early then begin
    let prod = ref 1.0 and iprod = ref 0 in
    for i = 1 to 2 * !mu do
      prod := !prod *. sqrt (!nu +. !dmu +. 1.0 -. float_of_int i);
      let r = dxadj env { x = !prod; ix = !iprod } in
      prod := r.x;
      iprod := r.ix
    done;
    for i = !j to l do
      let c1 = !prod *. sqrt (!nu +. 0.5) in
      let v = pqa.(i - 1) in
      pqa.(i - 1) <- dxadj env { x = v.x *. c1; ix = v.ix + !iprod };
      if nu2 -. nu1 > 0.5 then begin
        prod := sqrt (!nu +. !dmu +. 1.0) *. !prod;
        if !nu <> !dmu -. 1.0 then prod := !prod /. sqrt (!nu -. !dmu +. 1.0);
        let r = dxadj env { x = !prod; ix = !iprod } in
        prod := r.x;
        iprod := r.ix;
        nu := !nu +. 1.0
      end
      else if !dmu >= !nu then begin
        prod := 0.0;
        iprod := 0;
        mu := !mu + 1;
        dmu := !dmu +. 1.0
      end
      else begin
        prod := sqrt (!nu +. !dmu +. 1.0) *. !prod;
        if !nu > !dmu then prod := !prod *. sqrt (!nu -. !dmu);
        let r = dxadj env { x = !prod; ix = !iprod } in
        prod := r.x;
        iprod := r.ix;
        mu := !mu + 1;
        dmu := !dmu +. 1.0
      end
    done
  end

let dxqmu env nu1 nu2 mu1 mu2 theta x sx kind pqa =
  dxpqnu env nu1 nu2 0 theta kind pqa;
  let pq2 = ref pqa.(0) in
  dxpqnu env nu1 nu2 1 theta kind pqa;
  let nu = nu1 in
  let k = ref 0 in
  let mu = ref 1 and dmu = ref 1.0 in
  let pq1 = ref pqa.(0) in
  let fin = ref false in
  if mu1 <= 0 then begin
    incr k;
    pqa.(!k - 1) <- !pq2;
    if mu2 < 1 then fin := true
  end;
  if (not !fin) && mu1 <= 1 then begin
    incr k;
    pqa.(!k - 1) <- !pq1;
    if mu2 <= 1 then fin := true
  end;
  if not !fin then begin
    let stop = ref false and it = ref 0 in
    let bound = mu2 + 2 in
    while (not !stop) && !it <= bound do
      incr it;
      let x1 = -2.0 *. !dmu *. x *. sx *. !pq1.x in
      let x2 = (nu +. !dmu) *. (nu -. !dmu +. 1.0) *. !pq2.x in
      let s = dxadd env { x = x1; ix = !pq1.ix } { x = -.x2; ix = !pq2.ix } in
      let p = dxadj env s in
      pq2 := !pq1;
      pq1 := p;
      mu := !mu + 1;
      dmu := !dmu +. 1.0;
      if !mu >= mu1 then begin
        incr k;
        pqa.(!k - 1) <- p;
        if mu2 <= !mu then stop := true
      end
    done;
    if not !stop then invalid_arg "dxqmu: mu recurrence hit its bound"
  end

(* DXQNU walks off both ends. PQA(K-1) at K = 1, and the backward loop misses
   NU1 on a fractional dnu1 and writes one too far, so both are bounded here. *)
let dxqnu env nu1 nu2 mu1 theta x sx kind pqa =
  let k = ref 0 in
  let pq2 = ref zero and pql2 = ref zero in
  let fin = ref false in
  if mu1 <> 1 then begin
    dxpqnu env nu1 nu2 0 theta kind pqa;
    if mu1 = 0 then fin := true
    else begin
      k := int_of_float (nu2 -. nu1 +. 1.5);
      pq2 := pqa.(!k - 1);
      if !k >= 2 then pql2 := pqa.(!k - 2)
    end
  end;
  if not !fin then begin
    dxpqnu env nu1 nu2 1 theta kind pqa;
    if mu1 <> 1 then begin
      let nu = ref nu2 in
      let pq1 = ref pqa.(!k - 1) in
      let pql1 = ref (if !k >= 2 then pqa.(!k - 2) else zero) in
      let dmu = ref 1.0 in
      let back = ref false and pass = ref 0 in
      while (not !fin) && (not !back) && !pass < 3 do
        incr pass;
        let mu = ref 1 in
        dmu := 1.0;
        let p = ref zero in
        let stop = ref false and it = ref 0 in
        while (not !stop) && !it <= mu1 + 2 do
          incr it;
          let x1 = -2.0 *. !dmu *. x *. sx *. !pq1.x in
          let x2 = (!nu +. !dmu) *. (!nu -. !dmu +. 1.0) *. !pq2.x in
          let s =
            dxadd env { x = x1; ix = !pq1.ix } { x = -.x2; ix = !pq2.ix }
          in
          p := dxadj env s;
          pq2 := !pq1;
          pq1 := !p;
          mu := !mu + 1;
          dmu := !dmu +. 1.0;
          if !mu >= mu1 then stop := true
        done;
        if not !stop then invalid_arg "dxqnu: mu recurrence hit its bound";
        pqa.(!k - 1) <- !p;
        if !k = 1 then fin := true
        else if !nu < nu2 then back := true
        else begin
          nu := !nu -. 1.0;
          pq2 := !pql2;
          pq1 := !pql1;
          k := !k - 1
        end
      done;
      if (not !fin) && not !back then
        invalid_arg "dxqnu: mu-wise passes hit their bound";
      if !back then begin
        let pq1 = ref pqa.(!k - 1) and pq2 = ref pqa.(!k) in
        let it = ref 0 in
        let bound = Array.length pqa + 2 in
        while !nu > nu1 +. 0.5 && !it <= bound do
          incr it;
          k := !k - 1;
          let x1 = (2.0 *. !nu +. 1.0) *. x *. !pq1.x /. (!nu +. !dmu) in
          let x2 = -.(!nu -. !dmu +. 1.0) *. !pq2.x /. (!nu +. !dmu) in
          let s = dxadd env { x = x1; ix = !pq1.ix } { x = x2; ix = !pq2.ix } in
          let p = dxadj env s in
          pq2 := !pq1;
          pq1 := p;
          pqa.(!k - 1) <- p;
          nu := !nu -. 1.0
        done;
        if !it > bound || !k <> 1 then
          invalid_arg "dxqnu: nu recurrence did not fill the vector"
      end
    end
  end

let dxlegf env dnu1 nudiff mu1 mu2 theta kind =
  let pi2 = 2.0 *. atan 1.0 in
  if nudiff < 0 then invalid_arg "dxlegf: nudiff < 0, argument 2";
  if dnu1 < -0.5 then invalid_arg "dxlegf: dnu1 < -0.5, argument 1";
  if mu2 < mu1 then invalid_arg "dxlegf: mu2 < mu1, argument 4";
  if mu1 < 0 then invalid_arg "dxlegf: mu1 < 0, argument 3";
  if theta <= 0.0 || theta > pi2 then
    invalid_arg "dxlegf: theta out of range, argument 5";
  if mu1 <> mu2 && nudiff > 0 then
    invalid_arg "dxlegf: nudiff > 0 and mu2 > mu1, argument 2";
  let l = mu2 - mu1 + nudiff + 1 in
  let pqa = Array.make l zero in
  let dnu2 = dnu1 +. float_of_int nudiff in
  let frac = Float.rem dnu1 1.0 <> 0.0 in
  if frac then begin
    match kind with
    | Pnorm ->
      invalid_arg "dxlegf: normalised P needs an integer dnu1, argument 1"
    | Pneg | Q | Ppos -> ()
  end;
  let early =
    if frac then false
    else
      match kind with
      | Ppos | Pnorm -> float_of_int mu1 > dnu2
      | Pneg | Q -> false
  in
  if early then pqa
  else begin
    let x = cos theta and sx = 1.0 /. sin theta in
    (match kind with
    | Q ->
      if mu2 = mu1 then dxqnu env dnu1 dnu2 mu1 theta x sx kind pqa
      else dxqmu env dnu1 dnu2 mu1 mu2 theta x sx kind pqa
    | Pneg | Ppos | Pnorm ->
      if mu2 - mu1 <= 0 then dxpqnu env dnu1 dnu2 mu1 theta kind pqa
      else dxpmu env dnu1 dnu2 mu1 mu2 theta x sx kind pqa;
      (match kind with
      | Ppos -> dxpmup env dnu1 dnu2 mu1 mu2 pqa
      | Pnorm -> dxpnrm env dnu1 dnu2 mu1 mu2 pqa
      | Pneg | Q -> ()));
    for i = 0 to l - 1 do
      pqa.(i) <- dxred env pqa.(i)
    done;
    pqa
  end

let dxnrmp env nu mu1 mu2 darg mode =
  if nu < 0 then invalid_arg "dxnrmp: nu < 0, argument 1";
  if mu1 < 0 then invalid_arg "dxnrmp: mu1 < 0, argument 2";
  if mu1 > mu2 then invalid_arg "dxnrmp: mu1 > mu2, argument 3";
  let k = mu2 - mu1 + 1 in
  let dpn = Array.make k zero in
  let special =
    if nu = 0 then true
    else
      match mode with
      | By_x ->
        if Float.abs darg > 1.0 then
          invalid_arg "dxnrmp: darg out of range, argument 4";
        Float.abs darg = 1.0
      | By_theta ->
        if Float.abs darg > 4.0 *. atan 1.0 then
          invalid_arg "dxnrmp: darg out of range, argument 4";
        darg = 0.0
  in
  if special then
    if mu1 > 0 then (dpn, 0)
    else begin
      let v = sqrt (float_of_int nu +. 0.5) in
      let negate =
        nu mod 2 <> 0
        && not (match mode with By_x -> darg = 1.0 | By_theta -> true)
      in
      dpn.(0) <- { x = (if negate then -.v else v); ix = 0 };
      (dpn, 1)
    end
  else begin
    let sx, tx, isig =
      match mode with
      | By_x ->
        let xa = Float.abs darg in
        let sx = sqrt ((1.0 +. xa) *. (0.5 -. xa +. 0.5)) in
        let tx = darg /. sx in
        ( sx,
          tx,
          int_of_float (log10 (2.0 *. float_of_int nu *. (5.0 +. (tx *. tx))))
        )
      | By_theta ->
        let sx = Float.abs (sin darg) in
        let tx = cos darg /. sx in
        ( sx,
          tx,
          int_of_float
            (log10 (2.0 *. float_of_int nu *. (5.0 +. Float.abs (darg *. tx))))
        )
    in
    let mu = ref mu2 and i = ref k in
    while !i > 0 && !mu > nu do
      dpn.(!i - 1) <- zero;
      decr i;
      decr mu
    done;
    if !i = 0 then (dpn, 0)
    else begin
      let mu = ref nu in
      let p1 = ref zero and p2 = ref { x = 1.0; ix = 0 } in
      let p3 = ref 0.5 and dk = ref 2.0 in
      for _j = 1 to nu do
        p3 := (!dk +. 1.0) /. !dk *. !p3;
        p2 := dxadj env { x = !p2.x *. sx; ix = !p2.ix };
        dk := !dk +. 2.0
      done;
      p2 := dxadj env { x = !p2.x *. sqrt !p3; ix = !p2.ix };
      let s = 2.0 *. tx and t = 1.0 /. float_of_int nu in
      let fin = ref false in
      if mu2 >= nu then begin
        dpn.(!i - 1) <- !p2;
        decr i;
        if !i = 0 then fin := true
      end;
      let it = ref 0 in
      let bound = nu + 2 in
      while (not !fin) && !it <= bound do
        incr it;
        let p = float_of_int !mu *. t in
        let c1 = 1.0 /. sqrt ((1.0 -. p +. t) *. (1.0 +. p)) in
        let c2 = s *. p *. c1 *. !p2.x in
        let c1 = -.sqrt ((1.0 +. p +. t) *. (1.0 -. p)) *. c1 *. !p1.x in
        let pv = dxadd env { x = c2; ix = !p2.ix } { x = c1; ix = !p1.ix } in
        mu := !mu - 1;
        if !mu <= mu2 then begin
          dpn.(!i - 1) <- pv;
          decr i;
          if !i = 0 then fin := true
        end;
        if not !fin then begin
          p1 := !p2;
          p2 := pv;
          if !mu <= mu1 then fin := true
        end
      done;
      if not !fin then invalid_arg "dxnrmp: mu recurrence hit its bound";
      for idx = 0 to k - 1 do
        dpn.(idx) <- dxred env dpn.(idx)
      done;
      (dpn, isig)
    end
  end
