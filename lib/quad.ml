(* Reeve, Copyright 2026 Zane Hambly.
   Quadrature, ported from the SLATEC double precision routines. Where the
   Fortran took an EXTERNAL integrand these take a [float -> float] closure,
   and a matrix of Taylor derivatives is a [float array] in column major
   order with an explicit leading dimension, so the element (i, j) zero based
   lives at i + j * ldc. *)

open Mach

(* ---- Machine constants ---- *)

let nbits =
  int_of_float (log10_radix_dp *. float_of_int digits_dp /. 0.30102000)

(* ---- Tabulated data ---- *)

let davint x y n xlo xup =
  let ans = ref 0.0 and ierr = ref 1 in
  let exception Done in
  (try
     if xlo > xup then begin
       ierr := 2;
       raise Done
     end;
     if xlo = xup then raise Done;
     if n < 2 then begin
       ierr := 5;
       raise Done
     end;
     let i = ref 2 and scanning = ref true in
     while !scanning && !i <= n do
       if x.(!i - 1) <= x.(!i - 2) then begin
         ierr := 4;
         raise Done
       end;
       if x.(!i - 1) > xup then scanning := false else i := !i + 1
     done;
     if n < 3 then begin
       let slope = (y.(1) -. y.(0)) /. (x.(1) -. x.(0)) in
       let fl = y.(0) +. (slope *. (xlo -. x.(0))) in
       let fr = y.(1) +. (slope *. (xup -. x.(1))) in
       ans := 0.5 *. (fl +. fr) *. (xup -. xlo);
       raise Done
     end;
     if x.(n - 3) < xlo then begin
       ierr := 3;
       raise Done
     end;
     if x.(2) > xup then begin
       ierr := 3;
       raise Done
     end;
     let i = ref 1 and guard = ref 0 in
     while x.(!i - 1) < xlo do
       incr guard;
       if !guard > n then
         failwith "davint: the search for the lower limit ran past the table";
       i := !i + 1
     done;
     let inlft = !i in
     let i = ref n and guard = ref 0 in
     while x.(!i - 1) > xup do
       incr guard;
       if !guard > n then
         failwith "davint: the search for the upper limit ran past the table";
       i := !i - 1
     done;
     let inrt = !i in
     if inrt - inlft < 2 then begin
       ierr := 3;
       raise Done
     end;
     let istart = if inlft = 1 then 2 else inlft in
     let istop = if inrt = n then n - 1 else inrt in
     let r3 = 3.0 and rp5 = 0.5 in
     let summ = ref 0.0 in
     let syl = ref xlo in
     let syl2 = ref (!syl *. !syl) in
     let syl3 = ref (!syl2 *. !syl) in
     let ca = ref 0.0 and cb = ref 0.0 and cc = ref 0.0 in
     for i = istart to istop do
       let x1 = x.(i - 2) and x2 = x.(i - 1) and x3 = x.(i) in
       let x12 = x1 -. x2 and x13 = x1 -. x3 and x23 = x2 -. x3 in
       let term1 = y.(i - 2) /. (x12 *. x13) in
       let term2 = -.(y.(i - 1) /. (x12 *. x23)) in
       let term3 = y.(i) /. (x13 *. x23) in
       let a = term1 +. term2 +. term3 in
       let b =
         (-.(x2 +. x3) *. term1)
         -. ((x1 +. x3) *. term2)
         -. ((x1 +. x2) *. term3)
       in
       let c =
         (x2 *. x3 *. term1) +. (x1 *. x3 *. term2) +. (x1 *. x2 *. term3)
       in
       if i > istart then begin
         ca := 0.5 *. (a +. !ca);
         cb := 0.5 *. (b +. !cb);
         cc := 0.5 *. (c +. !cc)
       end
       else begin
         ca := a;
         cb := b;
         cc := c
       end;
       let syu = x2 in
       let syu2 = syu *. syu in
       let syu3 = syu2 *. syu in
       summ :=
         !summ
         +. (!ca *. (syu3 -. !syl3) /. r3)
         +. (!cb *. rp5 *. (syu2 -. !syl2))
         +. (!cc *. (syu -. !syl));
       ca := a;
       cb := b;
       cc := c;
       syl := syu;
       syl2 := syu2;
       syl3 := syu3
     done;
     let syu = xup in
     ans :=
       !summ
       +. (!ca *. ((syu *. syu *. syu) -. !syl3) /. r3)
       +. (!cb *. rp5 *. ((syu *. syu) -. !syl2))
       +. (!cc *. (syu -. !syl))
   with Done -> ());
  (!ans, !ierr)

(* ---- Breakpoint search and the piecewise polynomial ---- *)

let dintrv xt lxt x ilo =
  let ilo = ref ilo and ileft = ref 1 and mflag = ref 0 in
  let exception Done in
  (try
     let ihi = ref (!ilo + 1) in
     if !ihi >= lxt then begin
       if x >= xt.(lxt - 1) then begin
         mflag := 1;
         ileft := lxt;
         raise Done
       end;
       if lxt <= 1 then begin
         mflag := -1;
         ileft := 1;
         raise Done
       end;
       ilo := lxt - 1;
       ihi := lxt
     end;
     if x >= xt.(!ihi - 1) then begin
       let istep = ref 1 and guard = ref 0 and searching = ref true in
       while !searching do
         incr guard;
         if !guard > lxt + 2 then
           failwith "dintrv: the upward search did not reach the table end";
         ilo := !ihi;
         ihi := !ilo + !istep;
         if !ihi >= lxt then begin
           if x >= xt.(lxt - 1) then begin
             mflag := 1;
             ileft := lxt;
             raise Done
           end;
           ihi := lxt;
           searching := false
         end
         else if x < xt.(!ihi - 1) then searching := false
         else istep := !istep * 2
       done
     end
     else begin
       if x >= xt.(!ilo - 1) then begin
         mflag := 0;
         ileft := !ilo;
         raise Done
       end;
       let istep = ref 1 and guard = ref 0 and searching = ref true in
       while !searching do
         incr guard;
         if !guard > lxt + 2 then
           failwith
             "dintrv: the downward search did not reach the table start";
         ihi := !ilo;
         ilo := !ihi - !istep;
         if !ilo <= 1 then begin
           ilo := 1;
           if x >= xt.(0) then searching := false
           else begin
             mflag := -1;
             ileft := 1;
             raise Done
           end
         end
         else if x >= xt.(!ilo - 1) then searching := false
         else istep := !istep * 2
       done
     end;
     let narrow = ref true and guard = ref 0 in
     while !narrow do
       incr guard;
       if !guard > lxt + 2 then
         failwith "dintrv: the bisection did not narrow the interval";
       let middle = (!ilo + !ihi) / 2 in
       if middle = !ilo then begin
         mflag := 0;
         ileft := !ilo;
         narrow := false
       end
       else if x < xt.(middle - 1) then ihi := middle
       else ilo := middle
     done
   with Done -> ());
  (!ilo, !ileft, !mflag)

let dppval ldc c xi lxi k ideriv x =
  if k < 1 then invalid_arg "dppval: k must be at least 1";
  if ldc < k then invalid_arg "dppval: ldc must be at least k";
  if lxi < 1 then invalid_arg "dppval: lxi must be at least 1";
  if ideriv < 0 || ideriv >= k then
    invalid_arg "dppval: ideriv must satisfy 0 <= ideriv < k";
  let _, i, _ = dintrv xi lxi x 1 in
  let kk = ref (k - ideriv) in
  let dx = x -. xi.(i - 1) in
  let v = ref 0.0 and j = ref k in
  while !kk > 0 do
    v := (!v /. float_of_int !kk *. dx) +. c.(!j - 1 + ((i - 1) * ldc));
    j := !j - 1;
    kk := !kk - 1
  done;
  !v

(* ---- The adaptive eight point Legendre-Gauss rule ---- *)

let gx1 = 1.83434642495649805e-01
let gx2 = 5.25532409916328986e-01
let gx3 = 7.96666477413626740e-01
let gx4 = 9.60289856497536232e-01
let gw1 = 3.62683783378361983e-01
let gw2 = 3.13706645877887287e-01
let gw3 = 2.22381034453374471e-01
let gw4 = 1.01228536290376259e-01
let g8_sq2 = 1.41421356
let g8_kmx = 5000
let g8_kml = 6
let g8_max_steps = 20_000

let g8 f x h =
  let s1 = f (x -. (gx1 *. h)) +. f (x +. (gx1 *. h)) in
  let s2 = f (x -. (gx2 *. h)) +. f (x +. (gx2 *. h)) in
  let s3 = f (x -. (gx3 *. h)) +. f (x +. (gx3 *. h)) in
  let s4 = f (x -. (gx4 *. h)) +. f (x +. (gx4 *. h)) in
  h *. (((gw1 *. s1) +. (gw2 *. s2)) +. ((gw3 *. s3) +. (gw4 *. s4)))

let gaus8 name kern a b err =
  let nlmx = min 60 (nbits * 5 / 8) in
  let aa = Array.make 61 0.0
  and gr = Array.make 61 0.0
  and hh = Array.make 61 0.0
  and vl = Array.make 61 0.0
  and lr = Array.make 61 0 in
  let ans = ref 0.0 and ierr = ref 1 and ce = ref 0.0 in
  let exception Done in
  (try
     if a = b then raise Done;
     let lmx = ref nlmx in
     if b <> 0.0 && Float.copy_sign 1.0 b *. a > 0.0 then begin
       let c = Float.abs (1.0 -. (a /. b)) in
       if c <= 0.1 then begin
         if c <= 0.0 then raise Done;
         let anib = 0.5 -. (log c /. 0.69314718) in
         let nib = int_of_float anib in
         lmx := min nlmx (nbits - nib - 7);
         if !lmx < 1 then begin
           ierr := -1;
           raise Done
         end
       end
     end;
     let tol =
       if err = 0.0 then sqrt eps_dp
       else
         Float.max (Float.abs err) (2.0 ** float_of_int (5 - nbits)) /. 2.0
     in
     hh.(1) <- (b -. a) /. 4.0;
     aa.(1) <- a;
     lr.(1) <- 1;
     let l = ref 1 in
     let est = ref (kern (aa.(1) +. (2.0 *. hh.(1))) (2.0 *. hh.(1))) in
     let k = ref 8 in
     let area = ref (Float.abs !est) in
     let ef = ref 0.5 in
     let mxl = ref 0 in
     let eps = ref tol in
     let glr = ref 0.0 and vr = ref 0.0 in
     let steps = ref 0 in
     while true do
       let refining = ref true in
       while !refining do
         incr steps;
         if !steps > g8_max_steps then
           failwith (name ^ ": the subdivision did not terminate");
         let gl = kern (aa.(!l) +. hh.(!l)) hh.(!l) in
         gr.(!l) <- kern (aa.(!l) +. (3.0 *. hh.(!l))) hh.(!l);
         k := !k + 16;
         area :=
           !area +. (Float.abs gl +. Float.abs gr.(!l) -. Float.abs !est);
         glr := gl +. gr.(!l);
         let ee = Float.abs (!est -. !glr) *. !ef in
         let ae = Float.max (!eps *. !area) (tol *. Float.abs !glr) in
         if ee <= ae then refining := false
         else begin
           if !k > g8_kmx then lmx := g8_kml;
           if !l >= !lmx then begin
             mxl := 1;
             refining := false
           end
           else begin
             l := !l + 1;
             eps := !eps *. 0.5;
             ef := !ef /. g8_sq2;
             hh.(!l) <- hh.(!l - 1) *. 0.5;
             lr.(!l) <- -1;
             aa.(!l) <- aa.(!l - 1);
             est := gl
           end
         end
       done;
       ce := !ce +. (!est -. !glr);
       if lr.(!l) <= 0 then vl.(!l) <- !glr
       else begin
         vr := !glr;
         let up = ref true in
         while !up && !l > 1 do
           l := !l - 1;
           eps := !eps *. 2.0;
           ef := !ef *. g8_sq2;
           if lr.(!l) <= 0 then begin
             vl.(!l) <- vl.(!l + 1) +. !vr;
             up := false
           end
           else vr := vl.(!l + 1) +. !vr
         done;
         if !up then begin
           ans := !vr;
           if !mxl <> 0 && Float.abs !ce > 2.0 *. tol *. !area then ierr := 2;
           raise Done
         end
       end;
       est := gr.(!l - 1);
       lr.(!l) <- 1;
       aa.(!l) <- aa.(!l) +. (4.0 *. hh.(!l))
     done
   with Done -> ());
  (!ans, !ierr, if err < 0.0 then !ce else err)

(* LMN was write-only in the Fortran, its one read commented out at
   dgaus8.f:144. dqnc79 still uses it. *)
let dgaus8 f a b err = gaus8 "dgaus8" (g8 f) a b err

let dppgq8 f ldc c xi lxi kk id a b err =
  if kk < 1 then invalid_arg "dppgq8: kk must be at least 1";
  if ldc < kk then invalid_arg "dppgq8: ldc must be at least kk";
  if lxi < 1 then invalid_arg "dppgq8: lxi must be at least 1";
  if id < 0 || id >= kk then
    invalid_arg "dppgq8: id must satisfy 0 <= id < kk";
  gaus8 "dppgq8" (g8 (fun t -> f t *. dppval ldc c xi lxi kk id t)) a b err

(* ---- The adaptive seven point Newton-Cotes rule ---- *)

type qnc_step = Refine | Deeper | Rightward

let qw1 = 41.0 /. 140.0
let qw2 = 216.0 /. 140.0
let qw3 = 27.0 /. 140.0
let qw4 = 272.0 /. 140.0
let qnc_sq2 = sqrt 2.0
let qnc_kml = 7
let qnc_kmx = 5000
let qnc_nlmn = 2
let qnc_max_steps = 50_000

(* K is undefined on the early returns in the Fortran, assigned at line 140 but
   branched past at 117. This hands back 0. *)
let dqnc79 f a b err =
  let nlmx = min 99 (nbits * 4 / 5) in
  let aa = Array.make 100 0.0
  and hh = Array.make 100 0.0
  and q7r = Array.make 100 0.0
  and vl = Array.make 100 0.0
  and f1 = Array.make 100 0.0
  and f2 = Array.make 100 0.0
  and f3 = Array.make 100 0.0
  and f4 = Array.make 100 0.0
  and f5 = Array.make 100 0.0
  and f6 = Array.make 100 0.0
  and f7 = Array.make 100 0.0
  and lr = Array.make 100 0 in
  let fv = Array.make 14 0.0 in
  let ans = ref 0.0 and ierr = ref 1 and ce = ref 0.0 and kount = ref 0 in
  let exception Done in
  (try
     if a = b then begin
       ierr := -1;
       raise Done
     end;
     let lmx = ref nlmx and lmn = ref qnc_nlmn in
     if b <> 0.0 && Float.copy_sign 1.0 b *. a > 0.0 then begin
       let c = Float.abs (1.0 -. (a /. b)) in
       if c <= 0.1 then begin
         if c <= 0.0 then begin
           ierr := -1;
           raise Done
         end;
         let nib = int_of_float (0.5 -. (log c /. log 2.0)) in
         lmx := min nlmx (nbits - nib - 4);
         if !lmx < 2 then begin
           ierr := -1;
           raise Done
         end;
         lmn := min !lmn !lmx
       end
     end;
     let tol =
       if err = 0.0 then sqrt eps_dp
       else Float.max (Float.abs err) (2.0 ** float_of_int (5 - nbits))
     in
     hh.(1) <- (b -. a) /. 12.0;
     aa.(1) <- a;
     lr.(1) <- 1;
     let i = ref 1 in
     while !i <= 11 do
       fv.(!i) <- f (a +. (float_of_int (!i - 1) *. hh.(1)));
       i := !i + 2
     done;
     let blocal = b in
     fv.(13) <- f blocal;
     kount := 7;
     let l = ref 1 in
     let area = ref 0.0 and q7 = ref 0.0 and ef = ref (256.0 /. 255.0) in
     let bank = ref 0.0 and eps = ref tol in
     let q7l = ref 0.0 and q13 = ref 0.0 and vr = ref 0.0 in
     let state = ref Refine in
     let steps = ref 0 in
     while true do
       incr steps;
       if !steps > qnc_max_steps then
         failwith "dqnc79: the subdivision did not terminate";
       match !state with
       | Refine ->
         let i = ref 2 in
         while !i <= 12 do
           fv.(!i) <- f (aa.(!l) +. (float_of_int (!i - 1) *. hh.(!l)));
           i := !i + 2
         done;
         kount := !kount + 6;
         q7l :=
           hh.(!l)
           *. (((qw1 *. (fv.(1) +. fv.(7))) +. (qw2 *. (fv.(2) +. fv.(6))))
              +. ((qw3 *. (fv.(3) +. fv.(5))) +. (qw4 *. fv.(4))));
         q7r.(!l) <-
           hh.(!l)
           *. (((qw1 *. (fv.(7) +. fv.(13))) +. (qw2 *. (fv.(8) +. fv.(12))))
              +. ((qw3 *. (fv.(9) +. fv.(11))) +. (qw4 *. fv.(10))));
         area :=
           !area +. (Float.abs !q7l +. Float.abs q7r.(!l) -. Float.abs !q7);
         if !l < !lmn then state := Deeper
         else begin
           q13 := !q7l +. q7r.(!l);
           let ee = Float.abs (!q7 -. !q13) *. !ef in
           let ae = !eps *. !area in
           let test = Float.min (ae +. (0.8 *. !bank)) (10.0 *. ae) in
           let test =
             Float.max
               (Float.max test (tol *. Float.abs !q13))
               (0.00003 *. tol *. !area)
           in
           let descend = ref false in
           if ee <= test then ce := !ce +. ((!q7 -. !q13) /. 255.0)
           else begin
             if !kount > qnc_kmx then lmx := min qnc_kml !lmx;
             if !l < !lmx then descend := true
             else ce := !ce +. (!q7 -. !q13)
           end;
           if !descend then state := Deeper
           else begin
             bank := !bank +. (ae -. ee);
             if !bank < 0.0 then bank := 0.0;
             if lr.(!l) <= 0 then begin
               vl.(!l) <- !q13;
               state := Rightward
             end
             else begin
               vr := !q13;
               let up = ref true in
               while !up && !l > 1 do
                 if !l <= 17 then ef := !ef *. qnc_sq2;
                 eps := !eps *. 2.0;
                 l := !l - 1;
                 if lr.(!l) <= 0 then begin
                   vl.(!l) <- vl.(!l + 1) +. !vr;
                   up := false;
                   state := Rightward
                 end
                 else vr := vl.(!l + 1) +. !vr
               done;
               if !up then begin
                 ans := !vr;
                 if Float.abs !ce > 2.0 *. tol *. !area then ierr := 2;
                 raise Done
               end
             end
           end
         end
       | Deeper ->
         l := !l + 1;
         eps := !eps *. 0.5;
         if !l <= 17 then ef := !ef /. qnc_sq2;
         hh.(!l) <- hh.(!l - 1) *. 0.5;
         lr.(!l) <- -1;
         aa.(!l) <- aa.(!l - 1);
         q7 := !q7l;
         f1.(!l) <- fv.(7);
         f2.(!l) <- fv.(8);
         f3.(!l) <- fv.(9);
         f4.(!l) <- fv.(10);
         f5.(!l) <- fv.(11);
         f6.(!l) <- fv.(12);
         f7.(!l) <- fv.(13);
         fv.(13) <- fv.(7);
         fv.(11) <- fv.(6);
         fv.(9) <- fv.(5);
         fv.(7) <- fv.(4);
         fv.(5) <- fv.(3);
         fv.(3) <- fv.(2);
         state := Refine
       | Rightward ->
         q7 := q7r.(!l - 1);
         lr.(!l) <- 1;
         aa.(!l) <- aa.(!l) +. (12.0 *. hh.(!l));
         fv.(1) <- f1.(!l);
         fv.(3) <- f2.(!l);
         fv.(5) <- f3.(!l);
         fv.(7) <- f4.(!l);
         fv.(9) <- f5.(!l);
         fv.(11) <- f6.(!l);
         fv.(13) <- f7.(!l);
         state := Refine
     done
   with Done -> ());
  (!ans, !ierr, !kount)

(* ---- Products with a piecewise polynomial ---- *)

let dpfqad f ldc c xi lxi k id x1 x2 tol =
  if k < 1 then invalid_arg "dpfqad: k must be at least 1";
  if ldc < k then invalid_arg "dpfqad: ldc must be at least k";
  if id < 0 || id >= k then invalid_arg "dpfqad: id must satisfy 0 <= id < k";
  if lxi < 1 then invalid_arg "dpfqad: lxi must be at least 1";
  let wtol = Float.max eps_dp 1.0e-18 in
  if not (tol >= wtol && tol <= 0.1) then
    invalid_arg "dpfqad: tol below the unit roundoff or above 0.1";
  let ierr = ref 1 in
  let aa = Float.min x1 x2 and bb = Float.max x1 x2 in
  if aa = bb then (0.0, !ierr)
  else begin
    let ilo, il1, _ = dintrv xi lxi aa 1 in
    let _, il2, _ = dintrv xi lxi bb ilo in
    let q = ref 0.0 in
    for left = il1 to il2 do
      let ta = xi.(left - 1) in
      let a = if left = 1 then aa else Float.max aa ta in
      let tb = if left < lxi then xi.(left) else bb in
      let b = Float.min bb tb in
      let ans, iflg, _ = dppgq8 f ldc c xi lxi k id a b tol in
      if iflg > 1 then ierr := 2;
      q := !q +. ans
    done;
    ((if x1 > x2 then -. !q else !q), !ierr)
  end

module Raw = struct
  let davint = davint
  let dgaus8 = dgaus8
  let dqnc79 = dqnc79
  let dppgq8 = dppgq8
  let dpfqad = dpfqad
end

(* ---- The OCaml surface ---- *)

type result = { value : float; error : float option; converged : bool }

exception Too_narrow of float * float

let default_tol = sqrt eps_dp

let signed_tol tol =
  let t = Float.abs tol in
  -.(if t = 0.0 then default_tol else t)

let gauss8 ?(tol = default_tol) f a b =
  let ans, ierr, est = dgaus8 f a b (signed_tol tol) in
  if ierr < 0 then raise (Too_narrow (a, b));
  { value = ans; error = Some (Float.abs est); converged = ierr = 1 }

let newton_cotes7 ?(tol = default_tol) f a b =
  if a = b then { value = 0.0; error = None; converged = true }
  else begin
    let ans, ierr, _ = dqnc79 f a b (signed_tol tol) in
    if ierr < 0 then raise (Too_narrow (a, b));
    { value = ans; error = None; converged = ierr = 1 }
  end

let integrate_table x y lo up =
  let n = Array.length x in
  if n < 2 then
    invalid_arg "Quad.integrate_table: fewer than two tabulated points";
  if Array.length y < n then
    invalid_arg "Quad.integrate_table: fewer values than abscissas";
  if up < lo then
    invalid_arg "Quad.integrate_table: the upper limit is below the lower";
  let ans, ierr = davint x y n lo up in
  if ierr = 3 then
    invalid_arg
      "Quad.integrate_table: fewer than three abscissas between the limits";
  if ierr = 4 then
    invalid_arg "Quad.integrate_table: the abscissas do not increase";
  if ierr <> 1 then invalid_arg "Quad.integrate_table: improper input";
  ans
