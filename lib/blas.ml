(* Reeve, Copyright 2026 Zane Hambly.
   BLAS level 1, ported from the reference Fortran. A vector is a
   [float array], a base offset and an increment, and a negative increment
   walks the vector backwards from its far end as it does in the Fortran.
   Those shapes are also re-exported as Raw, with an OCaml surface over whole
   arrays under it. *)

open Mach

(* ---- Strides ---- *)

let start n inc = if inc < 0 then (1 - n) * inc else 0

(* ---- Level 1 ---- *)

let daxpy n da dx xo incx dy yo incy =
  if n > 0 && da <> 0.0 then
    if incx = 1 && incy = 1 then
      for i = 0 to n - 1 do
        dy.(yo + i) <- dy.(yo + i) +. (da *. dx.(xo + i))
      done
    else begin
      let ix = start n incx and iy = start n incy in
      for i = 0 to n - 1 do
        let jy = yo + iy + (i * incy) in
        dy.(jy) <- dy.(jy) +. (da *. dx.(xo + ix + (i * incx)))
      done
    end

let dscal n da dx xo incx =
  if n > 0 && incx > 0 && da <> 1.0 then
    for i = 0 to n - 1 do
      dx.(xo + (i * incx)) <- da *. dx.(xo + (i * incx))
    done

let dcopy n dx xo incx dy yo incy =
  if n > 0 then begin
    let ix = start n incx and iy = start n incy in
    for i = 0 to n - 1 do
      dy.(yo + iy + (i * incy)) <- dx.(xo + ix + (i * incx))
    done
  end

let dswap n dx xo incx dy yo incy =
  if n > 0 then begin
    let ix = start n incx and iy = start n incy in
    for i = 0 to n - 1 do
      let jx = xo + ix + (i * incx) and jy = yo + iy + (i * incy) in
      let t = dx.(jx) in
      dx.(jx) <- dy.(jy);
      dy.(jy) <- t
    done
  end

let drot n dx xo incx dy yo incy dc ds =
  if n > 0 then begin
    let ix = start n incx and iy = start n incy in
    for i = 0 to n - 1 do
      let jx = xo + ix + (i * incx) and jy = yo + iy + (i * incy) in
      let t = (dc *. dx.(jx)) +. (ds *. dy.(jy)) in
      dy.(jy) <- (dc *. dy.(jy)) -. (ds *. dx.(jx));
      dx.(jx) <- t
    done
  end

let drotg da db =
  let ada = Float.abs da and adb = Float.abs db in
  let roe = if ada > adb then da else db in
  let scale = ada +. adb in
  if scale = 0.0 then (0.0, 0.0, 1.0, 0.0)
  else begin
    let u = da /. scale and v = db /. scale in
    let r = Float.copy_sign (scale *. sqrt ((u *. u) +. (v *. v))) roe in
    let c = da /. r and s = db /. r in
    let z =
      if ada > adb then s
      else if c <> 0.0 then 1.0 /. c
      else 1.0
    in
    (r, z, c, s)
  end

let ddot n dx xo incx dy yo incy =
  if n <= 0 then 0.0
  else if incx = 1 && incy = 1 then begin
    let m = n mod 5 in
    let s = ref 0.0 in
    for i = 0 to m - 1 do
      s := !s +. (dx.(xo + i) *. dy.(yo + i))
    done;
    if n < 5 then !s
    else begin
      let i = ref m in
      while !i < n do
        s :=
          !s
          +. (dx.(xo + !i) *. dy.(yo + !i))
          +. (dx.(xo + !i + 1) *. dy.(yo + !i + 1))
          +. (dx.(xo + !i + 2) *. dy.(yo + !i + 2))
          +. (dx.(xo + !i + 3) *. dy.(yo + !i + 3))
          +. (dx.(xo + !i + 4) *. dy.(yo + !i + 4));
        i := !i + 5
      done;
      !s
    end
  end
  else begin
    let ix = start n incx and iy = start n incy in
    let s = ref 0.0 in
    for i = 0 to n - 1 do
      s := !s +. (dx.(xo + ix + (i * incx)) *. dy.(yo + iy + (i * incy)))
    done;
    !s
  end

let dasum n dx xo incx =
  if n <= 0 || incx <= 0 then 0.0
  else if incx = 1 then begin
    let m = n mod 6 in
    let s = ref 0.0 in
    for i = 0 to m - 1 do
      s := !s +. Float.abs dx.(xo + i)
    done;
    if n < 6 then !s
    else begin
      let i = ref m in
      while !i < n do
        s :=
          !s
          +. Float.abs dx.(xo + !i)
          +. Float.abs dx.(xo + !i + 1)
          +. Float.abs dx.(xo + !i + 2)
          +. Float.abs dx.(xo + !i + 3)
          +. Float.abs dx.(xo + !i + 4)
          +. Float.abs dx.(xo + !i + 5);
        i := !i + 6
      done;
      !s
    end
  end
  else begin
    let s = ref 0.0 in
    for i = 0 to n - 1 do
      s := !s +. Float.abs dx.(xo + (i * incx))
    done;
    !s
  end

(* ---- The scaled sum of squares behind dnrm2 ---- *)

let dnrm2 n dx xo incx =
  if n <= 0 then 0.0
  else begin
    let notbig = ref true in
    let asml = ref 0.0 and amed = ref 0.0 and abig = ref 0.0 in
    let ix = ref (start n incx) in
    for _i = 1 to n do
      let ax = Float.abs dx.(xo + !ix) in
      if ax > tbig then begin
        let t = ax *. sbig in
        abig := !abig +. (t *. t);
        notbig := false
      end
      else if ax < tsml then begin
        if !notbig then begin
          let t = ax *. ssml in
          asml := !asml +. (t *. t)
        end
      end
      else amed := !amed +. (ax *. ax);
      ix := !ix + incx
    done;
    let scl, sumsq =
      if !abig > 0.0 then begin
        if !amed > 0.0 || !amed <> !amed then
          abig := !abig +. (!amed *. sbig *. sbig);
        (1.0 /. sbig, !abig)
      end
      else if !asml > 0.0 then
        if !amed > 0.0 || !amed <> !amed then begin
          let a = sqrt !amed and s = sqrt !asml /. ssml in
          let ymin, ymax = if s > a then (a, s) else (s, a) in
          let r = ymin /. ymax in
          (1.0, ymax *. ymax *. (1.0 +. (r *. r)))
        end
        else (1.0 /. ssml, !asml)
      else (1.0, !amed)
    in
    scl *. sqrt sumsq
  end

let idamax n dx xo incx =
  if n < 1 || incx <= 0 then -1
  else begin
    let l = ref 0 and v = ref (Float.abs dx.(xo)) in
    for i = 1 to n - 1 do
      let w = Float.abs dx.(xo + (i * incx)) in
      if w > !v then begin
        v := w;
        l := i
      end
    done;
    !l
  end

module Raw = struct
  let daxpy = daxpy
  let dscal = dscal
  let dcopy = dcopy
  let dswap = dswap
  let drot = drot
  let drotg = drotg
  let ddot = ddot
  let dasum = dasum
  let dnrm2 = dnrm2
  let idamax = idamax
end

(* ---- The OCaml surface ---- *)

let count name inc len =
  if inc < 1 then
    invalid_arg (name ^ ": the increment has to be positive");
  (len + inc - 1) / inc

let pair name incx x incy y =
  let nx = count name incx (Array.length x) in
  let ny = count name incy (Array.length y) in
  if nx <> ny then
    invalid_arg (name ^ ": the vectors disagree on the number of elements");
  nx

let axpy ?(incx = 1) ?(incy = 1) alpha x y =
  daxpy (pair "Blas.axpy" incx x incy y) alpha x 0 incx y 0 incy

let dot ?(incx = 1) ?(incy = 1) x y =
  ddot (pair "Blas.dot" incx x incy y) x 0 incx y 0 incy

let nrm2 ?(inc = 1) x = dnrm2 (count "Blas.nrm2" inc (Array.length x)) x 0 inc
let asum ?(inc = 1) x = dasum (count "Blas.asum" inc (Array.length x)) x 0 inc

let scal ?(inc = 1) alpha x =
  dscal (count "Blas.scal" inc (Array.length x)) alpha x 0 inc

let copy ?(incx = 1) ?(incy = 1) x y =
  dcopy (pair "Blas.copy" incx x incy y) x 0 incx y 0 incy

let swap ?(incx = 1) ?(incy = 1) x y =
  dswap (pair "Blas.swap" incx x incy y) x 0 incx y 0 incy

let rot ?(incx = 1) ?(incy = 1) x y c s =
  drot (pair "Blas.rot" incx x incy y) x 0 incx y 0 incy c s

let rotg = drotg

let iamax ?(inc = 1) x =
  let n = count "Blas.iamax" inc (Array.length x) in
  if n < 1 then invalid_arg "Blas.iamax: the vector is empty";
  idamax n x 0 inc
