let knu_x1 = 2.0
let knu_x2 = 17.0
let knu_pi = 3.14159265358979e+00
let knu_rthpi = 1.25331413731550
let knu_miller_limit = 160

let dbsknu x fnu kode n y =
  let elim = 2.303 *. (float_of_int (-min_exp_dp) *. log10_radix_dp -. 3.0) in
  let tol = Float.max eps_2_dp 1.0e-15 in
  if x <= 0.0 then invalid_arg "dbsknu: x <= 0";
  if fnu < 0.0 then invalid_arg "dbsknu: order fnu < 0";
  if kode < 1 || kode > 2 then invalid_arg "dbsknu: kode not 1 or 2";
  if n < 1 then invalid_arg "dbsknu: n < 1";
  let rx = 2.0 /. x in
  let inu = ref (int_of_float (fnu +. 0.5)) in
  let dnu = fnu -. float_of_int !inu in
  let half_odd = Float.abs dnu = 0.5 in
  let dnu2 = if half_odd || Float.abs dnu < tol then 0.0 else dnu *. dnu in
  let series () =
    let a1 = 1.0 -. dnu and a2 = 1.0 +. dnu in
    let t1 = 1.0 /. dgamma_pos a1 and t2 = 1.0 /. dgamma_pos a2 in
    let g1 =
      if Float.abs dnu > 0.1 then (t1 -. t2) /. (dnu +. dnu)
      else begin
        let s = ref cc_tab.(0) in
        let ak = ref 1.0 in
        let k = ref 2 and brk = ref false in
        while not !brk && !k <= 8 do
          ak := !ak *. dnu2;
          let tm = cc_tab.(!k - 1) *. !ak in
          s := !s +. tm;
          if Float.abs tm < tol then brk := true;
          incr k
        done;
        -. !s
      end
    in
    let g2 = (t1 +. t2) *. 0.5 in
    let flrx = log rx in
    let fmu = dnu *. flrx in
    let smu = ref 1.0 in
    let fc = ref 1.0 in
    if dnu <> 0.0 then begin
      fc := dnu *. knu_pi;
      fc := !fc /. sin !fc;
      if fmu <> 0.0 then smu := sinh fmu /. fmu
    end;
    let f = ref (!fc *. (g1 *. cosh fmu +. g2 *. flrx *. !smu)) in
    let fce = exp fmu in
    let p = ref (0.5 *. fce /. t2) in
    let q = ref (0.5 /. (fce *. t1)) in
    let ak = ref 1.0 in
    let ck = ref 1.0 in
    let bk = ref 1.0 in
    let s1 = ref !f in
    let s2 = ref !p in
    if !inu > 0 || n > 1 then begin
      if x >= tol then begin
        let cx = x *. x *. 0.25 in
        let brk = ref false and it = ref 0 in
        while not !brk do
          incr it;
          if !it > iter_limit then failwith "dbsknu: power series did not converge";
          f := (!ak *. !f +. !p +. !q) /. (!bk -. dnu2);
          p := !p /. (!ak -. dnu);
          q := !q /. (!ak +. dnu);
          ck := !ck *. cx /. !ak;
          let t1 = !ck *. !f in
          s1 := !s1 +. t1;
          let t2 = !ck *. (!p -. !ak *. !f) in
          s2 := !s2 +. t2;
          bk := !bk +. !ak +. !ak +. 1.0;
          ak := !ak +. 1.0;
          let s =
            Float.abs t1 /. (1.0 +. Float.abs !s1)
            +. Float.abs t2 /. (1.0 +. Float.abs !s2)
          in
          if s <= tol then brk := true
        done
      end;
      s2 := !s2 *. rx;
      if kode <> 1 then begin
        let e = exp x in
        s1 := !s1 *. e;
        s2 := !s2 *. e
      end;
      `Two (!s1, !s2)
    end
    else begin
      if x >= tol then begin
        let cx = x *. x *. 0.25 in
        let brk = ref false and it = ref 0 in
        while not !brk do
          incr it;
          if !it > iter_limit then failwith "dbsknu: power series did not converge";
          f := (!ak *. !f +. !p +. !q) /. (!bk -. dnu2);
          p := !p /. (!ak -. dnu);
          q := !q /. (!ak +. dnu);
          ck := !ck *. cx /. !ak;
          let t1 = !ck *. !f in
          s1 := !s1 +. t1;
          bk := !bk +. !ak +. !ak +. 1.0;
          ak := !ak +. 1.0;
          if Float.abs t1 /. (1.0 +. Float.abs !s1) <= tol then brk := true
        done
      end;
      `One (if kode = 1 then !s1 else !s1 *. exp x)
    end
  in
  let asymptotic_or_miller () =
    let iflag = ref 0 in
    let coef = ref (knu_rthpi /. sqrt x) in
    if kode <> 2 then begin
      if x > elim then iflag := 1 else coef := !coef *. exp (-.x)
    end;
    let s1, s2, skip =
      if half_odd then (!coef, !coef, false)
      else if x > knu_x2 then begin
        let nn = if !inu = 0 && n = 1 then 1 else 2 in
        let dnu2 = dnu +. dnu in
        let fmu = ref (if Float.abs dnu2 >= tol then dnu2 *. dnu2 else 0.0) in
        let ex = x *. 8.0 in
        let s1 = ref 0.0 and s2 = ref 0.0 in
        for _k = 1 to nn do
          s1 := !s2;
          let s = ref 1.0 in
          let ak = ref 0.0 in
          let ck = ref 1.0 in
          let sqk = ref 1.0 in
          let dk = ref ex in
          let j = ref 1 and brk = ref false in
          while not !brk && !j <= 30 do
            ck := !ck *. (!fmu -. !sqk) /. !dk;
            s := !s +. !ck;
            dk := !dk +. ex;
            ak := !ak +. 8.0;
            sqk := !sqk +. !ak;
            if Float.abs !ck < tol then brk := true;
            incr j
          done;
          s2 := !s *. !coef;
          fmu := !fmu +. 8.0 *. dnu +. 4.0
        done;
        if nn <= 1 then (!s2, !s2, true) else (!s1, !s2, false)
      end
      else begin
        let etest = cos (knu_pi *. dnu) /. (knu_pi *. x *. tol) in
        let a = Array.make knu_miller_limit 0.0 in
        let b = Array.make knu_miller_limit 0.0 in
        let fks = ref 1.0 in
        let fhs = ref 0.25 in
        let fk = ref 0.0 in
        let ck = ref (x +. x +. 2.0) in
        let p1 = ref 0.0 and p2 = ref 1.0 in
        let k = ref 0 in
        let result = ref (0.0, 0.0, false) in
        let found = ref false in
        while not !found do
          incr k;
          if !k > knu_miller_limit then
            failwith "dbsknu: miller algorithm did not converge";
          fk := !fk +. 1.0;
          let ak = (!fhs -. dnu2) /. (!fks +. !fk) in
          let bk = !ck /. (!fk +. 1.0) in
          let pt = !p2 in
          p2 := bk *. !p2 -. ak *. !p1;
          p1 := pt;
          a.(!k - 1) <- ak;
          b.(!k - 1) <- bk;
          ck := !ck +. 2.0;
          fks := !fks +. !fk +. !fk +. 1.0;
          fhs := !fhs +. !fk +. !fk;
          if etest <= !fk *. !p1 then begin
            let kk = ref !k in
            let s = ref 1.0 in
            p1 := 0.0;
            p2 := 1.0;
            for _i = 1 to !k do
              let pt = !p2 in
              p2 := (b.(!kk - 1) *. !p2 -. !p1) /. a.(!kk - 1);
              p1 := pt;
              s := !s +. !p2;
              decr kk
            done;
            let s1 = !coef *. (!p2 /. !s) in
            if !inu <= 0 && n <= 1 then result := (s1, s1, true)
            else result := (s1, s1 *. (x +. dnu +. 0.5 -. !p1 /. !p2) /. x, false);
            found := true
          end
        done;
        !result
      end
    in
    (s1, s2, skip, !iflag)
  in
  let start = if (not half_odd) && x <= knu_x1 then series () else `Other in
  match start with
  | `One v ->
    y.(0) <- v;
    0
  | _ ->
    let s1, s2, skip, iflag =
      match start with
      | `Two (a, b) -> (a, b, false, 0)
      | _ -> asymptotic_or_miller ()
    in
    let s1 = ref s1 and s2 = ref s2 in
    let ck = ref 0.0 in
    if not skip then begin
      ck := (dnu +. dnu +. 2.0) /. x;
      if n = 1 then decr inu;
      if !inu > 0 then begin
        for _i = 1 to !inu do
          let st = !s2 in
          s2 := !ck *. !s2 +. !s1;
          s1 := st;
          ck := !ck +. rx
        done;
        if n = 1 then s1 := !s2
      end
      else if n <= 1 then s1 := !s2
    end;
    if iflag <> 1 then begin
      y.(0) <- !s1;
      if n > 1 then begin
        y.(1) <- !s2;
        for i = 3 to n do
          y.(i - 1) <- !ck *. y.(i - 2) +. y.(i - 3);
          ck := !ck +. rx
        done
      end;
      0
    end
    else begin
      let nz = ref 1 in
      y.(0) <- 0.0;
      let s = -.x +. log !s1 in
      if s >= -.elim then begin
        y.(0) <- exp s;
        nz := 0
      end;
      if n = 1 then !nz
      else begin
        y.(1) <- 0.0;
        nz := !nz + 1;
        let s = -.x +. log !s2 in
        if s >= -.elim then begin
          nz := !nz - 1;
          y.(1) <- exp s
        end;
        if n = 2 then !nz
        else begin
          let kk = ref 2 in
          let tail = ref true in
          if !nz >= 2 then begin
            let i = ref 3 and came_on_scale = ref false in
            while (not !came_on_scale) && !i <= n do
              kk := !i;
              let st = !s2 in
              s2 := !ck *. !s2 +. !s1;
              s1 := st;
              ck := !ck +. rx;
              let s = -.x +. log !s2 in
              nz := !nz + 1;
              y.(!i - 1) <- 0.0;
              if s >= -.elim then begin
                y.(!i - 1) <- exp s;
                nz := !nz - 1;
                came_on_scale := true
              end;
              incr i
            done;
            if not !came_on_scale then tail := false
          end;
          if not !tail then !nz
          else if !kk = n then !nz
          else begin
            s2 := !s2 *. !ck +. !s1;
            ck := !ck +. rx;
            incr kk;
            y.(!kk - 1) <- exp (-.x +. log !s2);
            if !kk = n then !nz
            else begin
              incr kk;
              for i = !kk to n do
                y.(i - 1) <- !ck *. y.(i - 2) +. y.(i - 3);
                ck := !ck +. rx
              done;
              !nz
            end
          end
        end
      end
    end

let besk_nulim = [| 35; 70 |]

let dbesk x fnu kode n y =
  let elim = 2.303 *. (float_of_int (-min_exp_dp) *. log10_radix_dp -. 3.0) in
  let xlim = tiny_dp *. 1.0e+3 in
  let overflow () =
    invalid_arg "dbesk: overflow, fnu or n too large or x too small"
  in
  if kode < 1 || kode > 2 then invalid_arg "dbesk: kode not 1 or 2";
  if fnu < 0.0 then invalid_arg "dbesk: order fnu < 0";
  if x <= 0.0 then invalid_arg "dbesk: x <= 0";
  if x < xlim then overflow ();
  if n < 1 then invalid_arg "dbesk: n < 1";
  let etx = float_of_int (kode - 1) in
  let nd = ref n in
  let nud = ref (int_of_float fnu) in
  let dnu = fnu -. float_of_int !nud in
  let gnu = ref fnu in
  let nn = ref (min 2 !nd) in
  let fn = fnu +. float_of_int n -. 1.0 in
  let fnn = fn in
  let cn = ref 0.0 in
  let rtz = ref 0.0 in
  let underflow_test () =
    let zn = x /. !gnu in
    rtz := sqrt (1.0 +. zn *. zn);
    let gln = log ((1.0 +. !rtz) /. zn) in
    let t = !rtz *. (1.0 -. etx) +. etx /. (zn +. !rtz) in
    cn := -. !gnu *. (t -. gln)
  in
  let underflow_loop () =
    let res = ref `Done and brk = ref false in
    while not !brk do
      incr nud;
      decr nd;
      if !nd = 0 then begin
        brk := true;
        res := `Done
      end
      else begin
        nn := min 2 !nd;
        gnu := !gnu +. 1.0;
        if fnn >= 2.0 && !nud >= besk_nulim.(!nn - 1) then begin
          underflow_test ();
          if !cn >= -.elim then begin
            brk := true;
            res := `Asymp
          end
        end
      end
    done;
    !res
  in
  let do_asymp () =
    dasyik x !gnu kode (-1.0) rtz cn !nn y 0;
    if !nn <> 1 then begin
      let trx = 2.0 /. x in
      let tm = ref ((!gnu +. !gnu +. 2.0) /. x) in
      for i = 3 to !nd do
        y.(i - 1) <- !tm *. y.(i - 2) +. y.(i - 3);
        tm := !tm +. trx
      done
    end
  in
  let do_overflow_test () =
    if fn > 1.0 && -.fn *. (log x -. 0.693) > elim then overflow ();
    if dnu = 0.0 then begin
      let j = ref !nud in
      let stop = ref false in
      if !j <> 1 then begin
        incr j;
        y.(!j - 1) <- (if kode = 2 then dbesk0e x else dbesk0 x);
        if !nd = 1 then stop := true else incr j
      end;
      if not !stop then y.(!j - 1) <- (if kode = 2 then dbesk1e x else dbesk1 x)
    end
    else ignore (dbsknu x fnu kode !nd y)
  in
  let do_dbsknu () =
    let pair =
      if dnu <> 0.0 then begin
        let nb = if !nud = 0 && !nd = 1 then 1 else 2 in
        let w = Array.make 2 0.0 in
        ignore (dbsknu x dnu kode nb w);
        if nb = 1 then `Single w.(0) else `Both (w.(0), w.(1))
      end
      else begin
        let s1 = if kode = 2 then dbesk0e x else dbesk0 x in
        if !nud = 0 && !nd = 1 then `Single s1
        else `Both (s1, if kode = 2 then dbesk1e x else dbesk1 x)
      end
    in
    let trx = 2.0 /. x in
    let tm = ref 0.0 in
    let s1 = ref 0.0 and s2 = ref 0.0 in
    (match pair with
     | `Single v -> s1 := v
     | `Both (a, b) ->
       s1 := a;
       s2 := b;
       tm := (dnu +. dnu +. 2.0) /. x;
       if !nd = 1 then decr nud;
       if !nud > 0 then begin
         for _i = 1 to !nud do
           let s = !s2 in
           s2 := !tm *. !s2 +. !s1;
           s1 := s;
           tm := !tm +. trx
         done;
         if !nd = 1 then s1 := !s2
       end
       else if !nd <= 1 then s1 := !s2);
    y.(0) <- !s1;
    if !nd <> 1 then begin
      y.(1) <- !s2;
      for i = 3 to !nd do
        y.(i - 1) <- !tm *. y.(i - 2) +. y.(i - 3);
        tm := !tm +. trx
      done
    end
  in
  let route =
    if fn < 2.0 then
      if kode = 2 || x <= elim then `Overflow_test else `Underflow_update
    else begin
      let zn = x /. fn in
      if zn = 0.0 then overflow ();
      let r = sqrt (1.0 +. zn *. zn) in
      rtz := r;
      let gln = log ((1.0 +. r) /. zn) in
      let t = r *. (1.0 -. etx) +. etx /. (zn +. r) in
      cn := -.fn *. (t -. gln);
      if !cn > elim then overflow ();
      if !nud < besk_nulim.(!nn - 1) then
        if kode = 2 || x <= elim then `Via_dbsknu else `Underflow_update
      else if !nn = 1 then `Asymp
      else begin
        underflow_test ();
        `Asymp
      end
    end
  in
  let rec run route =
    match route with
    | `Underflow_update -> (
      match underflow_loop () with `Done -> () | `Asymp -> run `Asymp)
    | `Asymp ->
      if !cn < -.elim then (
        match underflow_loop () with `Done -> () | `Asymp -> do_asymp ())
      else do_asymp ()
    | `Overflow_test -> do_overflow_test ()
    | `Via_dbsknu -> do_dbsknu ()
  in
  run route;
  let nz = n - !nd in
  if nz <> 0 then begin
    for i = 1 to !nd do
      y.(n - i) <- y.(!nd - i)
    done;
    for i = 1 to nz do
      y.(i - 1) <- 0.0
    done
  end;
  nz
