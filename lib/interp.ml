(* Reeve, Copyright 2026 Zane Hambly.
   Polynomial interpolation, ported from the SLATEC double precision routines.
   [dplint] builds the Newton divided difference table, [dpolvl] evaluates it
   along with derivatives and [dpolcf] turns it into Taylor coefficients. *)

let dplint n x y c =
  if n <= 0 then invalid_arg "dplint: n is zero or negative";
  c.(0) <- y.(0);
  if n <> 1 then
    for k = 2 to n do
      c.(k - 1) <- y.(k - 1);
      let km1 = k - 1 in
      for i = 1 to km1 do
        let dif = x.(i - 1) -. x.(k - 1) in
        if dif = 0.0 then
          invalid_arg "dplint: the abscissas are not distinct";
        c.(k - 1) <- (c.(i - 1) -. c.(k - 1)) /. dif
      done
    done

let dpolcf xx n x c d work =
  for k = 1 to n do
    d.(k - 1) <- c.(k - 1)
  done;
  if n <> 1 then begin
    work.(0) <- 1.0;
    let pone = ref c.(0) and ptwo = ref c.(0) in
    let nm1 = n - 1 in
    for k = 2 to n do
      let km1 = k - 1 in
      let npkm1 = n + k - 1 in
      work.(npkm1 - 1) <- xx -. x.(km1 - 1);
      work.(k - 1) <- work.(npkm1 - 1) *. work.(km1 - 1);
      ptwo := !pone +. (work.(k - 1) *. c.(k - 1));
      pone := !ptwo
    done;
    d.(0) <- !ptwo;
    if n <> 2 then
      for k = 2 to nm1 do
        let km1 = k - 1 in
        let km2n = k - 2 + n in
        let nmkp1 = n - k + 1 in
        for i = 2 to nmkp1 do
          let km2npi = km2n + i in
          let im1 = i - 1 in
          let km1pi = km1 + i in
          work.(i - 1) <- (work.(km2npi - 1) *. work.(im1 - 1)) +. work.(i - 1);
          d.(k - 1) <- d.(k - 1) +. (work.(i - 1) *. d.(km1pi - 1))
        done
      done
  end

let dpolvl nder xx yp n x c work =
  let ierr = 1 in
  if nder <= 0 then begin
    let pione = ref 1.0 and pone = ref c.(0) in
    if n = 1 then (!pone, ierr)
    else begin
      let ptwo = ref c.(0) in
      for k = 2 to n do
        let pitwo = (xx -. x.(k - 2)) *. !pione in
        pione := pitwo;
        ptwo := !pone +. (pitwo *. c.(k - 1));
        pone := !ptwo
      done;
      (!ptwo, ierr)
    end
  end
  else if n > 1 then begin
    let izero, ndr = if nder < n then (0, nder) else (1, n - 1) in
    let m = ndr + 1 in
    let mm = ref m in
    for k = 1 to ndr do
      yp.(k - 1) <- c.(k)
    done;
    work.(0) <- 1.0;
    let pone = ref c.(0) and ptwo = ref c.(0) in
    for k = 2 to n do
      let km1 = k - 1 in
      let npkm1 = n + k - 1 in
      work.(npkm1 - 1) <- xx -. x.(km1 - 1);
      work.(k - 1) <- work.(npkm1 - 1) *. work.(km1 - 1);
      ptwo := !pone +. (work.(k - 1) *. c.(k - 1));
      pone := !ptwo
    done;
    let yfit = !ptwo in
    if n <> 2 then begin
      if m = n then mm := ndr;
      for k = 2 to !mm do
        let nmkp1 = n - k + 1 in
        let km1 = k - 1 in
        let km2pn = k - 2 + n in
        for i = 2 to nmkp1 do
          let km2pni = km2pn + i in
          let im1 = i - 1 in
          let km1pi = km1 + i in
          work.(i - 1) <- (work.(km2pni - 1) *. work.(im1 - 1)) +. work.(i - 1);
          yp.(km1 - 1) <- yp.(km1 - 1) +. (work.(i - 1) *. c.(km1pi - 1))
        done
      done;
      if ndr <> 1 then begin
        let fac = ref 1.0 in
        for k = 2 to ndr do
          let xk = float_of_int k in
          fac := xk *. !fac;
          yp.(k - 1) <- !fac *. yp.(k - 1)
        done
      end
    end;
    if izero <> 0 then
      for k = n to nder do
        yp.(k - 1) <- 0.0
      done;
    (yfit, ierr)
  end
  else begin
    for k = 1 to nder do
      yp.(k - 1) <- 0.0
    done;
    (c.(0), ierr)
  end

module Raw = struct
  let dplint = dplint
  let dpolvl = dpolvl
  let dpolcf = dpolcf
end

(* ---- The OCaml surface ---- *)

type t = { x : float array; c : float array }

let make x y =
  let n = Array.length x in
  if n = 0 then invalid_arg "Interp.make: no points";
  if Array.length y < n then
    invalid_arg "Interp.make: fewer values than abscissas";
  let c = Array.make n 0.0 in
  dplint n x y c;
  { x = Array.sub x 0 n; c }

let points t = Array.length t.x

let eval t xx =
  let yfit, _ = dpolvl 0 xx [||] (Array.length t.x) t.x t.c [||] in
  yfit

let derivatives t xx nder =
  if nder < 0 then invalid_arg "Interp.derivatives: negative order";
  if nder = 0 then (eval t xx, [||])
  else begin
    let n = Array.length t.x in
    let yp = Array.make nder 0.0 in
    let work = Array.make (2 * n) 0.0 in
    let yfit, _ = dpolvl nder xx yp n t.x t.c work in
    (yfit, yp)
  end

let coefficients t xx =
  let n = Array.length t.x in
  let d = Array.make n 0.0 in
  let work = Array.make (2 * n) 0.0 in
  dpolcf xx n t.x t.c d work;
  d
