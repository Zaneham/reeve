let besi_rttpi = 3.98942280401433e-01
let besi_inlim = 80

let dbesi x alpha kode n y =
  let tol = Float.max eps_2_dp 1.0e-15 in
  let elim = 2.303 *. (float_of_int (-min_exp_dp) *. log10_radix_dp -. 3.0) in
  let tolln =
    Float.min (2.303 *. log10_radix_dp *. float_of_int (digits_dp + 1)) 34.5388
  in
  if n < 1 then invalid_arg "dbesi: n < 1";
  if kode < 1 || kode > 2 then invalid_arg "dbesi: kode not 1 or 2";
  if x < 0.0 then invalid_arg "dbesi: x < 0";
  if alpha < 0.0 then invalid_arg "dbesi: order alpha < 0";
  if x = 0.0 then begin
    if alpha = 0.0 then begin
      y.(0) <- 1.0;
      for i = 2 to n do
        y.(i - 1) <- 0.0
      done
    end
    else
      for i = 1 to n do
        y.(i - 1) <- 0.0
      done;
    0
  end
  else begin
    let kt = ref (if n = 1 then 2 else 1) in
    let nn = ref n in
    let nz = ref 0 in
    let ialp = int_of_float alpha in
    let fni = ref (float_of_int (ialp + n - 1)) in
    let fnf = alpha -. float_of_int ialp in
    let dfn = ref (!fni +. fnf) in
    let fnu = !dfn in
    let in_ = ref 0 in
    let xo2 = x *. 0.5 in
    let sxo2 = xo2 *. xo2 in
    let etx = float_of_int (kode - 1) in
    let sx = etx *. x in
    let temp = Array.make 3 0.0 in
    let fn = ref 0.0 and fnp1 = ref 0.0 and xo2l = ref 0.0 in
    let is = ref !kt in
    let i1 = ref 0 in
    let km = ref 0 and ns = ref 0 in
    let arg = ref 0.0 and earg = ref 0.0 in
    let gln = ref 0.0 and ra = ref 0.0 in
    let etx8 = ref 0.0 in
    let steps = ref 0 in
    let tick () =
      incr steps;
      if !steps > iter_limit then failwith "dbesi: recursion did not terminate"
    in
    let rec st_underflow_test () =
      tick ();
      let z = x /. !fn in
      ra := sqrt (1.0 +. z *. z);
      gln := log ((1.0 +. !ra) /. z);
      let t = !ra *. (1.0 -. etx) +. etx /. (z +. !ra) in
      arg := !fn *. (t -. !gln);
      st_check_arg ()
    and st_check_arg () =
      tick ();
      if !arg >= -.elim then st_asymp_call ()
      else begin
        y.(!nn - 1) <- 0.0;
        decr nn;
        fni := !fni -. 1.0;
        dfn := !fni +. fnf;
        fn := !dfn;
        if !nn < 1 then st_return_nz ()
        else begin
          if !nn = 1 then begin
            kt := 2;
            is := 2
          end;
          st_underflow_test ()
        end
      end
    and st_update_params () =
      tick ();
      is := 2;
      fni := !fni -. 1.0;
      dfn := !fni +. fnf;
      fn := !dfn;
      if !i1 = 2 then begin
        nz := n - !nn;
        st_backward_recur ()
      end
      else begin
        let z = x /. !fn in
        ra := sqrt (1.0 +. z *. z);
        gln := log ((1.0 +. !ra) /. z);
        let t = !ra *. (1.0 -. etx) +. etx /. (z +. !ra) in
        arg := !fn *. (t -. !gln);
        st_asymp_call ()
      end
    and st_asymp_call () =
      tick ();
      i1 := max (abs (3 - !is)) 1;
      dasyik x !fn kode 1.0 ra arg !i1 temp (!is - 1);
      if !is = 1 then st_update_params ()
      else if !is = 2 then begin
        nz := n - !nn;
        st_backward_recur ()
      end
      else begin
        let t = 1.0 /. (!fn *. !ra) in
        let ain = tolln /. (!gln +. sqrt ((!gln *. !gln) +. (t *. tolln))) +. 1.5 in
        in_ := int_of_float ain;
        if !in_ <= besi_inlim then st_backward_norm ()
        else if !km <> 0 then begin
          temp.(0) <- temp.(2);
          in_ := !ns;
          kt := 1;
          i1 := 0;
          st_update_params ()
        end
        else begin
          y.(0) <- temp.(2);
          !nz
        end
      end
    and st_series_init () =
      tick ();
      gln := dgamln !fnp1;
      arg := (!fn *. !xo2l) -. !gln -. sx;
      if !arg < -.elim then st_series_underflow ()
      else begin
        earg := exp !arg;
        st_series_loop ()
      end
    and st_series_loop () =
      tick ();
      let s = ref 1.0 in
      if x >= tol then begin
        let ak = ref 3.0 and t2 = ref 1.0 and t = ref 1.0 and s1 = ref !fn in
        let k = ref 1 and brk = ref false in
        while (not !brk) && !k <= 17 do
          let s2 = !t2 +. !s1 in
          t := !t *. sxo2 /. s2;
          s := !s +. !t;
          if Float.abs !t < tol then brk := true
          else begin
            t2 := !t2 +. !ak;
            ak := !ak +. 2.0;
            s1 := !s1 +. !fn
          end;
          incr k
        done
      end;
      temp.(!is - 1) <- !s *. !earg;
      if !is = 2 then begin
        nz := n - !nn;
        st_backward_recur ()
      end
      else if !is = 3 then begin
        km := int_of_float (Float.max (3.0 -. !fn) 0.0);
        let tfn = !fn +. float_of_int !km in
        let t = (!gln +. tfn -. 0.9189385332 -. (0.0833333333 /. tfn)) /. (tfn +. 0.5) in
        let ta = !xo2l -. t in
        let tb = -.(1.0 -. (1.0 /. tfn)) /. tfn in
        let ain = tolln /. (-.ta +. sqrt ((ta *. ta) -. (tolln *. tb))) +. 1.5 in
        in_ := int_of_float ain + !km;
        st_backward_norm ()
      end
      else begin
        earg := !earg *. !fn /. xo2;
        fni := !fni -. 1.0;
        dfn := !fni +. fnf;
        fn := !dfn;
        is := 2;
        st_series_loop ()
      end
    and st_series_underflow () =
      tick ();
      y.(!nn - 1) <- 0.0;
      decr nn;
      fnp1 := !fn;
      fni := !fni -. 1.0;
      dfn := !fni +. fnf;
      fn := !dfn;
      if !nn < 1 then st_return_nz ()
      else begin
        if !nn = 1 then begin
          kt := 2;
          is := 2
        end;
        if sxo2 > !fnp1 then st_underflow_test ()
        else begin
          arg := !arg -. !xo2l +. log !fnp1;
          if !arg >= -.elim then st_series_init () else st_series_underflow ()
        end
      end
    and st_return_nz () =
      nz := n - !nn;
      !nz
    and st_backward_recur () =
      tick ();
      if !kt = 2 then begin
        y.(0) <- temp.(1);
        !nz
      end
      else begin
        let s1 = ref temp.(0) and s2 = ref temp.(1) in
        let trx = 2.0 /. x in
        let dtm = ref !fni in
        let tm = ref ((!dtm +. fnf) *. trx) in
        let finish () =
          let k = ref (!nn + 1) in
          for _i = 3 to !nn do
            decr k;
            y.(!k - 3) <- (!tm *. y.(!k - 2)) +. y.(!k - 1);
            dtm := !dtm -. 1.0;
            tm := (!dtm +. fnf) *. trx
          done;
          !nz
        in
        if !in_ = 0 then begin
          y.(!nn - 1) <- !s1;
          y.(!nn - 2) <- !s2;
          if !nn = 2 then !nz else finish ()
        end
        else begin
          for _i = 1 to !in_ do
            let s = !s2 in
            s2 := (!tm *. !s2) +. !s1;
            s1 := s;
            dtm := !dtm -. 1.0;
            tm := (!dtm +. fnf) *. trx
          done;
          y.(!nn - 1) <- !s1;
          if !nn = 1 then !nz
          else begin
            y.(!nn - 2) <- !s2;
            if !nn = 2 then !nz else finish ()
          end
        end
      end
    and st_asymp_x_init () =
      tick ();
      etx8 := 8.0 *. x;
      is := !kt;
      in_ := 0;
      fn := fnu;
      st_asymp_x_loop ()
    and st_asymp_x_loop () =
      tick ();
      let dx0 = !fni +. !fni in
      let tm =
        if !fni <> 0.0 || Float.abs fnf >= tol then 4.0 *. fnf *. (!fni +. !fni +. fnf)
        else 0.0
      in
      let dtm = dx0 *. dx0 in
      let s1 = ref !etx8 in
      let dx = ref (-.(dtm -. 1.0 +. tm) /. !etx8) in
      let t = ref !dx in
      let s = ref (1.0 +. !dx) in
      let atol = tol *. Float.abs !s in
      let s2 = ref 1.0 in
      let ak = ref 8.0 in
      let k = ref 1 and brk = ref false in
      while (not !brk) && !k <= 25 do
        s1 := !s1 +. !etx8;
        s2 := !s2 +. !ak;
        dx := dtm -. !s2;
        let ap = !dx +. tm in
        t := -. !t *. ap /. !s1;
        s := !s +. !t;
        if Float.abs !t <= atol then brk := true else ak := !ak +. 8.0;
        incr k
      done;
      temp.(!is - 1) <- !s *. !earg;
      if !is = 2 then st_backward_recur ()
      else begin
        is := 2;
        fni := !fni -. 1.0;
        dfn := !fni +. fnf;
        fn := !dfn;
        st_asymp_x_loop ()
      end
    and st_backward_norm () =
      tick ();
      let trx = 2.0 /. x in
      let dtm = ref (!fni +. float_of_int !in_) in
      let tm = ref ((!dtm +. fnf) *. trx) in
      let ta = ref 0.0 and tb = ref tol in
      let kk = ref 1 and brk = ref false in
      while not !brk do
        for _i = 1 to !in_ do
          let s = !tb in
          tb := (!tm *. !tb) +. !ta;
          ta := s;
          dtm := !dtm -. 1.0;
          tm := (!dtm +. fnf) *. trx
        done;
        if !kk <> 1 then brk := true
        else begin
          ta := (!ta /. !tb) *. temp.(2);
          tb := temp.(2);
          kk := 2;
          in_ := !ns;
          if !ns = 0 then brk := true
        end
      done;
      y.(!nn - 1) <- !tb;
      nz := n - !nn;
      if !nn = 1 then !nz
      else begin
        tb := (!tm *. !tb) +. !ta;
        let k = ref (!nn - 1) in
        y.(!k - 1) <- !tb;
        if !nn = 2 then !nz
        else begin
          dtm := !dtm -. 1.0;
          tm := (!dtm +. fnf) *. trx;
          let last = !k - 1 in
          for _i = 1 to last do
            y.(!k - 2) <- (!tm *. y.(!k - 1)) +. y.(!k);
            dtm := !dtm -. 1.0;
            tm := (!dtm +. fnf) *. trx;
            decr k
          done;
          !nz
        end
      end
    in
    if sxo2 <= fnu +. 1.0 then begin
      fn := fnu;
      fnp1 := !fn +. 1.0;
      xo2l := log xo2;
      is := !kt;
      if x <= 0.5 then st_series_init ()
      else begin
        ns := 0;
        fni := !fni +. float_of_int !ns;
        dfn := !fni +. fnf;
        fn := !dfn;
        fnp1 := !fn +. 1.0;
        is := !kt;
        if n - 1 + !ns > 0 then is := 3;
        st_series_init ()
      end
    end
    else if x <= 12.0 then begin
      xo2l := log xo2;
      ns := int_of_float (sxo2 -. fnu);
      fni := !fni +. float_of_int !ns;
      dfn := !fni +. fnf;
      fn := !dfn;
      fnp1 := !fn +. 1.0;
      is := !kt;
      if n - 1 + !ns > 0 then is := 3;
      st_series_init ()
    end
    else begin
      let fnt = Float.max 17.0 (0.55 *. fnu *. fnu) in
      if x >= fnt then begin
        earg := besi_rttpi /. sqrt x;
        if kode = 2 then st_asymp_x_init ()
        else if x > elim then invalid_arg "dbesi: overflow, x too large for kode = 1"
        else begin
          earg := !earg *. exp x;
          st_asymp_x_init ()
        end
      end
      else begin
        ns := int_of_float (Float.max (36.0 -. fnu) 0.0);
        fni := !fni +. float_of_int !ns;
        dfn := !fni +. fnf;
        fn := !dfn;
        is := !kt;
        km := n - 1 + !ns;
        if !km > 0 then is := 3;
        if kode = 2 then st_underflow_test ()
        else if alpha < 1.0 then
          if x <= elim then st_underflow_test ()
          else invalid_arg "dbesi: overflow, x too large for kode = 1"
        else begin
          let z = x /. alpha in
          ra := sqrt (1.0 +. z *. z);
          gln := log ((1.0 +. !ra) /. z);
          let t = !ra *. (1.0 -. etx) +. etx /. (z +. !ra) in
          arg := alpha *. (t -. !gln);
          if !arg > elim then invalid_arg "dbesi: overflow, x too large for kode = 1"
          else if !km <> 0 then st_underflow_test ()
          else st_check_arg ()
        end
      end
    end
  end

let besj_rtwo = 1.34839972492648e+00
let besj_pdf = 7.85398163397448e-01
let besj_rttp = 7.97884560802865e-01
let besj_pidt = 1.57079632679490
let besj_inlim = 150
let besj_fnulim = [| 100.0; 60.0 |]

let dbesj x alpha n y =
  let tol = Float.max eps_2_dp 1.0e-15 in
  let elim1 = -2.303 *. ((float_of_int min_exp_dp *. log10_radix_dp) +. 3.0) in
  let rtol = 1.0 /. tol in
  let slim = tiny_dp *. rtol *. 1.0e+3 in
  let tolln =
    Float.min (2.303 *. log10_radix_dp *. float_of_int (digits_dp + 1)) 34.5388
  in
  if n < 1 then invalid_arg "dbesj: n < 1";
  if x < 0.0 then invalid_arg "dbesj: x < 0";
  if alpha < 0.0 then invalid_arg "dbesj: order alpha < 0";
  if x = 0.0 then begin
    if alpha = 0.0 then begin
      y.(0) <- 1.0;
      for i = 2 to n do
        y.(i - 1) <- 0.0
      done
    end
    else
      for i = 1 to n do
        y.(i - 1) <- 0.0
      done;
    0
  end
  else begin
    let kt = ref (if n = 1 then 2 else 1) in
    let nn = ref n in
    let nz = ref 0 in
    let ialp = int_of_float alpha in
    let fni = ref (float_of_int (ialp + n - 1)) in
    let fnf = alpha -. float_of_int ialp in
    let dfn = ref (!fni +. fnf) in
    let fnu = !dfn in
    let xo2 = x *. 0.5 in
    let sxo2 = xo2 *. xo2 in
    let temp = Array.make 3 0.0 in
    let wk = Array.make 7 0.0 in
    let fn = ref 0.0 and fnp1 = ref 0.0 and xo2l = ref 0.0 in
    let is = ref !kt in
    let i1 = ref 0 in
    let in_ = ref 0 and km = ref 0 and ns = ref 0 in
    let arg = ref 0.0 and earg = ref 0.0 and gln = ref 0.0 in
    let coef = ref 0.0 and etx = ref 0.0 in
    let fidal = ref 0.0 and dalpha = ref 0.0 in
    let sa = ref 0.0 and sb = ref 0.0 in
    let ta = ref 0.0 and tb = ref 0.0 and ak = ref 0.0 in
    let dtm = ref 0.0 and tm = ref 0.0 and trx = ref 0.0 in
    let kk = ref 0 in
    let steps = ref 0 in
    let tick () =
      incr steps;
      if !steps > iter_limit then failwith "dbesj: recursion did not terminate"
    in
    let rec st_uniform () =
      tick ();
      i1 := max (abs (3 - !is)) 1;
      let iflw = dasyjy djairy x !fn 1.0 !i1 temp (!is - 1) wk in
      if iflw <> 0 then begin
        y.(!nn - 1) <- 0.0;
        decr nn;
        fni := !fni -. 1.0;
        dfn := !fni +. fnf;
        fn := !dfn;
        if !nn < 1 then st_return_nz ()
        else begin
          if !nn = 1 then begin
            kt := 2;
            is := 2
          end;
          st_uniform ()
        end
      end
      else if !is = 1 then begin
        is := 2;
        fni := !fni -. 1.0;
        dfn := !fni +. fnf;
        fn := !dfn;
        if !i1 <> 2 then st_uniform () else st_backward_recur ()
      end
      else if !is = 2 then st_backward_recur ()
      else begin
        gln := wk.(2) +. wk.(1);
        let tau =
          if wk.(5) > 30.0 then begin
            let t = 0.5 *. tolln /. wk.(3) in
            ((0.0493827160 *. t -. 0.1111111111) *. t +. 0.6666666667) *. t *. wk.(5)
          end
          else begin
            let rden = ((pp_tab.(3) *. wk.(5)) +. pp_tab.(2)) *. wk.(5) +. 1.0 in
            let rzden = pp_tab.(0) +. (pp_tab.(1) *. wk.(5)) in
            rzden /. rden
          end
        in
        let tbl =
          if wk.(0) < 0.10 then
            (1.259921049 +. ((0.1679894730 +. (0.0887944358 *. wk.(0))) *. wk.(0)))
            /. wk.(6)
          else !gln /. wk.(4)
        in
        in_ := int_of_float ((tau /. tbl) +. 1.5);
        if !in_ <= besj_inlim then st_backward_setup ()
        else begin
          temp.(0) <- temp.(2);
          kt := 1;
          is := 2;
          fni := !fni -. 1.0;
          dfn := !fni +. fnf;
          fn := !dfn;
          if !i1 <> 2 then st_uniform () else st_backward_recur ()
        end
      end
    and st_series () =
      tick ();
      gln := dgamln !fnp1;
      arg := (!fn *. !xo2l) -. !gln;
      if !arg < -.elim1 then st_underflow ()
      else begin
        earg := exp !arg;
        st_series_inner ()
      end
    and st_series_inner () =
      tick ();
      let s = ref 1.0 in
      if x >= tol then begin
        let akl = ref 3.0 and t2 = ref 1.0 and t = ref 1.0 and s1 = ref !fn in
        let k = ref 1 and brk = ref false in
        while (not !brk) && !k <= 17 do
          let s2 = !t2 +. !s1 in
          t := -. !t *. sxo2 /. s2;
          s := !s +. !t;
          if Float.abs !t < tol then brk := true
          else begin
            t2 := !t2 +. !akl;
            akl := !akl +. 2.0;
            s1 := !s1 +. !fn
          end;
          incr k
        done
      end;
      temp.(!is - 1) <- !s *. !earg;
      if !is = 2 then st_backward_recur ()
      else if !is = 3 then begin
        km := int_of_float (Float.max (3.0 -. !fn) 0.0);
        let tfn = !fn +. float_of_int !km in
        let t = (!gln +. tfn -. 0.9189385332 -. (0.0833333333 /. tfn)) /. (tfn +. 0.5) in
        let tta = !xo2l -. t in
        let ttb = -.(1.0 -. (1.5 /. tfn)) /. tfn in
        let akm = tolln /. (-.tta +. sqrt ((tta *. tta) -. (tolln *. ttb))) +. 1.5 in
        in_ := !km + int_of_float akm;
        st_backward_setup ()
      end
      else begin
        earg := !earg *. !fn /. xo2;
        fni := !fni -. 1.0;
        dfn := !fni +. fnf;
        fn := !dfn;
        is := 2;
        st_series_inner ()
      end
    and st_underflow () =
      tick ();
      y.(!nn - 1) <- 0.0;
      decr nn;
      fnp1 := !fn;
      fni := !fni -. 1.0;
      dfn := !fni +. fnf;
      fn := !dfn;
      if !nn < 1 then st_return_nz ()
      else begin
        if !nn = 1 then begin
          kt := 2;
          is := 2
        end;
        if sxo2 > !fnp1 then st_uniform ()
        else begin
          arg := !arg -. !xo2l +. log !fnp1;
          if !arg >= -.elim1 then begin
            fnp1 := !fn +. 1.0;
            gln := dgamln !fnp1;
            arg := (!fn *. !xo2l) -. !gln;
            earg := exp !arg;
            st_series_inner ()
          end
          else st_underflow ()
        end
      end
    and st_return_nz () =
      nz := n - !nn;
      !nz
    and st_backward_recur () =
      tick ();
      let early = ref false in
      if !ns = 0 then begin
        nz := n - !nn;
        if !kt = 2 then begin
          y.(0) <- temp.(1);
          early := true
        end
        else begin
          y.(!nn - 1) <- temp.(0);
          y.(!nn - 2) <- temp.(1);
          if !nn = 2 then early := true
        end
      end;
      if !early then !nz
      else begin
        trx := 2.0 /. x;
        dtm := !fni;
        tm := (!dtm +. fnf) *. !trx;
        ak := 1.0;
        ta := temp.(0);
        tb := temp.(1);
        if Float.abs !ta <= slim then begin
          ta := !ta *. rtol;
          tb := !tb *. rtol;
          ak := tol
        end;
        kk := 2;
        in_ := !ns - 1;
        if !in_ = 0 then st_backward_indexed ()
        else if !ns <> 0 then st_backward_unindexed ()
        else begin
          let k = ref (!nn - 2) in
          for _i = 3 to !nn do
            let s = !tb in
            tb := (!tm *. !tb) -. !ta;
            ta := s;
            y.(!k - 1) <- !tb *. !ak;
            dtm := !dtm -. 1.0;
            tm := (!dtm +. fnf) *. !trx;
            decr k
          done;
          !nz
        end
      end
    and st_asymp_x () =
      tick ();
      let going = ref true in
      while !going do
        let d = !fidal +. !fidal in
        let dtm0 = d *. d in
        let tm0 =
          if !fidal <> 0.0 || Float.abs fnf >= tol then
            4.0 *. fnf *. (!fidal +. !fidal +. fnf)
          else 0.0
        in
        let t2 = ref ((dtm0 -. 1.0 +. tm0) /. !etx) in
        let s2 = ref !t2 in
        let relb = tol *. Float.abs !t2 in
        let t1 = ref !etx in
        let s1 = ref 1.0 in
        let fl = ref 1.0 in
        let akl = ref 8.0 in
        let k = ref 1 and brk = ref false in
        while (not !brk) && !k <= 13 do
          t1 := !t1 +. !etx;
          fl := !fl +. !akl;
          let ap = dtm0 -. !fl +. tm0 in
          t2 := -. !t2 *. ap /. !t1;
          s1 := !s1 +. !t2;
          t1 := !t1 +. !etx;
          akl := !akl +. 8.0;
          fl := !fl +. !akl;
          let ap = dtm0 -. !fl +. tm0 in
          t2 := !t2 *. ap /. !t1;
          s2 := !s2 +. !t2;
          if Float.abs !t2 <= relb then brk := true else akl := !akl +. 8.0;
          incr k
        done;
        temp.(!is - 1) <- !coef *. ((!s1 *. !sb) -. (!s2 *. !sa));
        if !is = 2 then going := false
        else begin
          fidal := !fidal +. 1.0;
          dalpha := !fidal +. fnf;
          is := 2;
          let t = !sa in
          sa := -. !sb;
          sb := t
        end
      done;
      if !kt = 2 then begin
        y.(0) <- temp.(1);
        !nz
      end
      else begin
        let s1 = ref temp.(0) and s2 = ref temp.(1) in
        let tx = 2.0 /. x in
        let tmf = ref (!dalpha *. tx) in
        let early = ref false in
        if !in_ <> 0 then begin
          for _i = 1 to !in_ do
            let s = !s2 in
            s2 := (!tmf *. !s2) -. !s1;
            tmf := !tmf +. tx;
            s1 := s
          done;
          if !nn = 1 then begin
            y.(0) <- !s2;
            early := true
          end
          else begin
            let s = !s2 in
            s2 := (!tmf *. !s2) -. !s1;
            tmf := !tmf +. tx;
            s1 := s
          end
        end;
        if !early then !nz
        else begin
          y.(0) <- !s1;
          y.(1) <- !s2;
          for i = 3 to !nn do
            y.(i - 1) <- (!tmf *. y.(i - 2)) -. y.(i - 3);
            tmf := !tmf +. tx
          done;
          !nz
        end
      end
    and st_backward_setup () =
      tick ();
      dtm := !fni +. float_of_int !in_;
      trx := 2.0 /. x;
      tm := (!dtm +. fnf) *. !trx;
      ta := 0.0;
      tb := tol;
      kk := 1;
      ak := 1.0;
      st_backward_unindexed ()
    and st_backward_unindexed () =
      tick ();
      let brk = ref false in
      while not !brk do
        for _i = 1 to !in_ do
          let s = !tb in
          tb := (!tm *. !tb) -. !ta;
          ta := s;
          dtm := !dtm -. 1.0;
          tm := (!dtm +. fnf) *. !trx
        done;
        if !kk <> 1 then brk := true
        else begin
          let s = temp.(2) in
          let scale = !ta /. !tb in
          ta := s;
          tb := s;
          if Float.abs s <= slim then begin
            ta := !ta *. rtol;
            tb := !tb *. rtol;
            ak := tol
          end;
          ta := !ta *. scale;
          kk := 2;
          in_ := !ns;
          if !ns = 0 then brk := true
        end
      done;
      st_backward_indexed ()
    and st_backward_indexed () =
      tick ();
      y.(!nn - 1) <- !tb *. !ak;
      nz := n - !nn;
      if !nn = 1 then !nz
      else begin
        let s = !tb in
        tb := (!tm *. !tb) -. !ta;
        ta := s;
        y.(!nn - 2) <- !tb *. !ak;
        if !nn = 2 then !nz
        else begin
          dtm := !dtm -. 1.0;
          tm := (!dtm +. fnf) *. !trx;
          let k = ref (!nn - 2) in
          for _i = 3 to !nn do
            let s = !tb in
            tb := (!tm *. !tb) -. !ta;
            ta := s;
            y.(!k - 1) <- !tb *. !ak;
            dtm := !dtm -. 1.0;
            tm := (!dtm +. fnf) *. !trx;
            decr k
          done;
          !nz
        end
      end
    in
    if sxo2 <= fnu +. 1.0 then begin
      fn := fnu;
      fnp1 := !fn +. 1.0;
      xo2l := log xo2;
      is := !kt;
      if x <= 0.50 then st_series ()
      else begin
        ns := 0;
        fni := !fni +. float_of_int !ns;
        dfn := !fni +. fnf;
        fn := !dfn;
        fnp1 := !fn +. 1.0;
        is := !kt;
        if n - 1 + !ns > 0 then is := 3;
        st_series ()
      end
    end
    else begin
      let tlim = Float.max 20.0 fnu in
      if x > tlim then begin
        let rtx = sqrt x in
        let tau = besj_rtwo *. rtx in
        if fnu <= tau +. besj_fnulim.(!kt - 1) then begin
          let inv = int_of_float (alpha -. tau +. 2.0) in
          let idalp =
            if inv <= 0 then begin
              in_ := 0;
              ialp
            end
            else begin
              in_ := inv;
              kt := 1;
              ialp - inv - 1
            end
          in
          is := !kt;
          fidal := float_of_int idalp;
          dalpha := !fidal +. fnf;
          let a = x -. (besj_pidt *. !dalpha) -. besj_pdf in
          sa := sin a;
          sb := cos a;
          coef := besj_rttp /. rtx;
          etx := 8.0 *. x;
          st_asymp_x ()
        end
        else begin
          fn := fnu;
          is := !kt;
          st_uniform ()
        end
      end
      else if x > 12.0 then begin
        ns := int_of_float (Float.max (36.0 -. fnu) 0.0);
        fni := !fni +. float_of_int !ns;
        dfn := !fni +. fnf;
        fn := !dfn;
        is := !kt;
        if n - 1 + !ns > 0 then is := 3;
        st_uniform ()
      end
      else begin
        xo2l := log xo2;
        ns := int_of_float (sxo2 -. fnu) + 1;
        fni := !fni +. float_of_int !ns;
        dfn := !fni +. fnf;
        fn := !dfn;
        fnp1 := !fn +. 1.0;
        is := !kt;
        if n - 1 + !ns > 0 then is := 3;
        st_series ()
      end
    end
  end
