(* Reeve, Copyright 2026 Zane Hambly.
   The LAPACK auxiliary routines, ported from the reference double precision
   Fortran. In the Fortran shapes a matrix is a [float array] in column major
   order with a base offset and an explicit leading dimension, so the element
   (i, j) zero based lives at ao + i + j * lda. The OCaml surface at the
   bottom takes a Mat.t instead and allocates its own workspace. *)

open Mach

type machine_parameter =
  | Eps
  | Safe_min
  | Base
  | Precision
  | Digits
  | Round
  | Min_exponent
  | Underflow
  | Max_exponent
  | Overflow

type uplo = Upper | Lower | Full
type norm = Max_abs | One_norm | Inf_norm | Frobenius

let xerbla srname info =
  raise
    (Invalid_argument
       (Printf.sprintf "%s: parameter number %d had an illegal value" srname
          info))

let fbase = 2.0
let feps = epsilon_float *. 0.5

let fsafmin =
  let small = 1.0 /. max_float in
  if small >= min_float then small *. (1.0 +. feps) else min_float

let dlamch cmach =
  match cmach with
  | Eps -> feps
  | Safe_min -> fsafmin
  | Base -> fbase
  | Precision -> feps *. fbase
  | Digits -> float_of_int fdigits
  | Round -> 1.0
  | Min_exponent -> float_of_int fminexp
  | Underflow -> min_float
  | Max_exponent -> float_of_int fmaxexp
  | Overflow -> max_float

let dlaisnan (din1 : float) (din2 : float) = din1 <> din2
let disnan din = dlaisnan din din

let dlapy2 x y =
  let x_is_nan = disnan x and y_is_nan = disnan y in
  if y_is_nan then y
  else if x_is_nan then x
  else begin
    let hugeval = dlamch Overflow in
    let xabs = Float.abs x and yabs = Float.abs y in
    let w = Float.max xabs yabs and z = Float.min xabs yabs in
    if z = 0.0 || w > hugeval then w
    else begin
      let r = z /. w in
      w *. sqrt (1.0 +. (r *. r))
    end
  end

let dlapy3 x y z =
  let hugeval = dlamch Overflow in
  let xabs = Float.abs x and yabs = Float.abs y and zabs = Float.abs z in
  let w = Float.max xabs (Float.max yabs zabs) in
  if w = 0.0 || w > hugeval then xabs +. yabs +. zabs
  else begin
    let rx = xabs /. w and ry = yabs /. w and rz = zabs /. w in
    w *. sqrt ((rx *. rx) +. (ry *. ry) +. (rz *. rz))
  end

let dlassq n x xo incx scale sumsq =
  if disnan scale || disnan sumsq then (scale, sumsq)
  else begin
    let scl = ref scale and ssq = ref sumsq in
    if !ssq = 0.0 then scl := 1.0;
    if !scl = 0.0 then begin
      scl := 1.0;
      ssq := 0.0
    end;
    if n <= 0 then (!scl, !ssq)
    else begin
      let notbig = ref true in
      let asml = ref 0.0 and amed = ref 0.0 and abig = ref 0.0 in
      let ix = ref (if incx < 0 then -(n - 1) * incx else 0) in
      for _i = 1 to n do
        let ax = Float.abs x.(xo + !ix) in
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
      if !ssq > 0.0 then begin
        let ax = !scl *. sqrt !ssq in
        if ax > tbig then
          if !scl > 1.0 then begin
            scl := !scl *. sbig;
            abig := !abig +. (!scl *. (!scl *. !ssq))
          end
          else abig := !abig +. (!scl *. (!scl *. (sbig *. (sbig *. !ssq))))
        else if ax < tsml then begin
          if !notbig then
            if !scl < 1.0 then begin
              scl := !scl *. ssml;
              asml := !asml +. (!scl *. (!scl *. !ssq))
            end
            else asml := !asml +. (!scl *. (!scl *. (ssml *. (ssml *. !ssq))))
        end
        else amed := !amed +. (!scl *. (!scl *. !ssq))
      end;
      if !abig > 0.0 then begin
        if !amed > 0.0 || disnan !amed then
          abig := !abig +. (!amed *. sbig *. sbig);
        (1.0 /. sbig, !abig)
      end
      else if !asml > 0.0 then
        if !amed > 0.0 || disnan !amed then begin
          let a = sqrt !amed and s = sqrt !asml /. ssml in
          let ymin, ymax = if s > a then (a, s) else (s, a) in
          let r = ymin /. ymax in
          (1.0, ymax *. ymax *. (1.0 +. (r *. r)))
        end
        else (1.0 /. ssml, !asml)
      else (1.0, !amed)
    end
  end

let dlaswp n a ao lda k1 k2 ipiv ipo incx =
  if incx <> 0 then begin
    let ix0, i1, i2, inc =
      if incx > 0 then (k1, k1, k2, 1)
      else (k1 + ((k1 - k2) * incx), k2, k1, -1)
    in
    let n32 = n / 32 * 32 in
    if n32 <> 0 then begin
      let j = ref 0 in
      while !j < n32 do
        let ix = ref ix0 and i = ref i1 in
        while (inc > 0 && !i <= i2) || (inc < 0 && !i >= i2) do
          let ip = ipiv.(ipo + !ix) in
          if ip <> !i then
            for k = !j to !j + 31 do
              let p = ao + !i + (k * lda) and q = ao + ip + (k * lda) in
              let temp = a.(p) in
              a.(p) <- a.(q);
              a.(q) <- temp
            done;
          ix := !ix + incx;
          i := !i + inc
        done;
        j := !j + 32
      done
    end;
    if n32 <> n then begin
      let ix = ref ix0 and i = ref i1 in
      while (inc > 0 && !i <= i2) || (inc < 0 && !i >= i2) do
        let ip = ipiv.(ipo + !ix) in
        if ip <> !i then
          for k = n32 to n - 1 do
            let p = ao + !i + (k * lda) and q = ao + ip + (k * lda) in
            let temp = a.(p) in
            a.(p) <- a.(q);
            a.(q) <- temp
          done;
        ix := !ix + incx;
        i := !i + inc
      done
    end
  end

let dlaset uplo m n alpha beta a ao lda =
  (match uplo with
  | Upper ->
    for j = 2 to n do
      for i = 1 to min (j - 1) m do
        a.(ao + i - 1 + ((j - 1) * lda)) <- alpha
      done
    done
  | Lower ->
    for j = 1 to min m n do
      for i = j + 1 to m do
        a.(ao + i - 1 + ((j - 1) * lda)) <- alpha
      done
    done
  | Full ->
    for j = 1 to n do
      for i = 1 to m do
        a.(ao + i - 1 + ((j - 1) * lda)) <- alpha
      done
    done);
  for i = 1 to min m n do
    a.(ao + i - 1 + ((i - 1) * lda)) <- beta
  done

let dlacpy uplo m n a ao lda b bo ldb =
  match uplo with
  | Upper ->
    for j = 1 to n do
      for i = 1 to min j m do
        b.(bo + i - 1 + ((j - 1) * ldb)) <- a.(ao + i - 1 + ((j - 1) * lda))
      done
    done
  | Lower ->
    for j = 1 to n do
      for i = j to m do
        b.(bo + i - 1 + ((j - 1) * ldb)) <- a.(ao + i - 1 + ((j - 1) * lda))
      done
    done
  | Full ->
    for j = 1 to n do
      for i = 1 to m do
        b.(bo + i - 1 + ((j - 1) * ldb)) <- a.(ao + i - 1 + ((j - 1) * lda))
      done
    done

let dlange norm m n a ao lda work =
  if min m n = 0 then 0.0
  else
    match norm with
    | Max_abs ->
      let value = ref 0.0 in
      for j = 1 to n do
        for i = 1 to m do
          let temp = Float.abs a.(ao + i - 1 + ((j - 1) * lda)) in
          if !value < temp || disnan temp then value := temp
        done
      done;
      !value
    | One_norm ->
      let value = ref 0.0 in
      for j = 1 to n do
        let sum = ref 0.0 in
        for i = 1 to m do
          sum := !sum +. Float.abs a.(ao + i - 1 + ((j - 1) * lda))
        done;
        if !value < !sum || disnan !sum then value := !sum
      done;
      !value
    | Inf_norm ->
      for i = 1 to m do
        work.(i - 1) <- 0.0
      done;
      for j = 1 to n do
        for i = 1 to m do
          work.(i - 1) <-
            work.(i - 1) +. Float.abs a.(ao + i - 1 + ((j - 1) * lda))
        done
      done;
      let value = ref 0.0 in
      for i = 1 to m do
        let temp = work.(i - 1) in
        if !value < temp || disnan temp then value := temp
      done;
      !value
    | Frobenius ->
      let scale = ref 0.0 and sum = ref 1.0 in
      for j = 1 to n do
        let s, q = dlassq m a (ao + ((j - 1) * lda)) 1 !scale !sum in
        scale := s;
        sum := q
      done;
      !scale *. sqrt !sum

let ieeeck ispec zero one =
  let posinf1 = one /. zero in
  let neginf1 = -.one /. zero in
  let negzro = one /. (neginf1 +. one) in
  let neginf2 = one /. negzro in
  let newzro = negzro +. zero in
  let posinf2 = one /. newzro in
  let neginf3 = neginf2 *. posinf2 in
  let posinf3 = posinf2 *. posinf2 in
  if
    not
      (not (posinf1 <= one)
      && not (neginf1 >= zero)
      && not (negzro <> zero)
      && not (neginf2 >= zero)
      && not (newzro <> zero)
      && not (posinf2 <= one)
      && not (neginf3 >= zero)
      && not (posinf3 <= one))
  then 0
  else if ispec = 0 then 1
  else begin
    let nan1 = posinf3 +. neginf3 in
    let nan2 = posinf3 /. neginf3 in
    let nan3 = posinf3 /. posinf3 in
    let nan4 = posinf3 *. zero in
    let nan5 = neginf3 *. negzro in
    let nan6 = nan5 *. zero in
    if
      not
        (nan1 = nan1 || nan2 = nan2 || nan3 = nan3 || nan4 = nan4
        || nan5 = nan5 || nan6 = nan6)
    then 1
    else 0
  end

let fpad s n =
  let l = String.length s in
  if l >= n then String.sub s 0 n else s ^ String.make (n - l) ' '

let fstreq a b =
  let n = max (String.length a) (String.length b) in
  String.equal (fpad a n) (fpad b n)

let fupper s n =
  let t = Bytes.of_string (fpad s n) in
  let ic = Char.code (Bytes.get t 0) in
  if ic >= 97 && ic <= 122 then begin
    Bytes.set t 0 (Char.chr (ic - 32));
    for i = 1 to min 5 (n - 1) do
      let c = Char.code (Bytes.get t i) in
      if c >= 97 && c <= 122 then Bytes.set t i (Char.chr (c - 32))
    done
  end;
  Bytes.to_string t

let iparmq ispec name _opts _n ilo ihi _lwork =
  let inmin = 12
  and inwin = 13
  and inibl = 14
  and ishfts = 15
  and iacc22 = 16
  and icost = 17 in
  let nmin = 75
  and k22min = 14
  and kacmin = 14
  and nibble = 14
  and knwswp = 500
  and rcost = 10 in
  let nh = ihi - ilo + 1 in
  let ns =
    if ispec = ishfts || ispec = inwin || ispec = iacc22 then begin
      let s = ref 2 in
      if nh >= 30 then s := 4;
      if nh >= 60 then s := 10;
      if nh >= 150 then
        s :=
          max 10
            (nh / int_of_float (Float.round (log (float_of_int nh) /. log 2.0)));
      if nh >= 590 then s := 64;
      if nh >= 3000 then s := 128;
      if nh >= 6000 then s := 256;
      max 2 (!s - (!s mod 2))
    end
    else 0
  in
  if ispec = inmin then nmin
  else if ispec = inibl then nibble
  else if ispec = ishfts then ns
  else if ispec = inwin then if nh <= knwswp then ns else 3 * ns / 2
  else if ispec = iacc22 then begin
    let subnam = fupper name 6 in
    let s26 = String.sub subnam 1 5 in
    if String.equal s26 "GGHRD" || String.equal s26 "GGHD3" then
      if nh >= k22min then 2 else 1
    else if String.equal (String.sub subnam 3 3) "EXC" then
      if nh >= k22min then 2 else if nh >= kacmin then 1 else 0
    else if
      String.equal s26 "HSEQR" || String.equal (String.sub subnam 1 4) "LAQR"
    then if ns >= k22min then 2 else if ns >= kacmin then 1 else 0
    else 0
  end
  else if ispec = icost then rcost
  else -1

(* The QP3RK branch never fires. ilaenv.f compares a four character slice
   against a five character literal and Fortran pads the short side, so
   DGEQP3RK falls through to the defaults. Kept that way. *)
let ilaenv ispec name opts n1 n2 n3 n4 =
  if ispec >= 1 && ispec <= 3 then begin
    let subnam = fupper name 16 in
    let c1 = subnam.[0] in
    let sname = Char.equal c1 'S' || Char.equal c1 'D' in
    let cname = Char.equal c1 'C' || Char.equal c1 'Z' in
    if not (cname || sname) then 1
    else begin
      let c2 = String.sub subnam 1 2 in
      let c3 = String.sub subnam 3 3 in
      let c4 = String.sub c3 1 2 in
      let twostage = Char.equal subnam.[10] '2' in
      let unitary = (sname && String.equal c2 "OR") || (cname && String.equal c2 "UN") in
      let gm () = Char.equal c3.[0] 'G' || Char.equal c3.[0] 'M' in
      let qrset () =
        match c4 with
        | "QR" | "RQ" | "LQ" | "QL" | "HR" | "TR" | "BR" -> true
        | _ -> false
      in
      let factored () =
        String.equal c3 "QRF" || String.equal c3 "RQF" || String.equal c3 "LQF"
        || String.equal c3 "QLF"
      in
      let qp3rk () = fstreq (String.sub subnam 3 4) "QP3RK" in
      let tall () = if n1 * n2 <= 131072 || n1 <= 8192 then n1 else 32768 / n2 in
      if ispec = 1 then begin
        let nb = ref 1 in
        if String.equal (String.sub subnam 1 5) "LAORH" then nb := 32
        else if String.equal c2 "GE" then begin
          if String.equal c3 "TRF" then nb := 64
          else if factored () then nb := 32
          else if String.equal c3 "QR " then nb := (if n3 = 1 then tall () else 1)
          else if String.equal c3 "LQ " then nb := (if n3 = 2 then tall () else 1)
          else if String.equal c3 "HRD" then nb := 32
          else if String.equal c3 "BRD" then nb := 32
          else if String.equal c3 "TRI" then nb := 64
          else if qp3rk () then nb := 32
        end
        else if String.equal c2 "PO" then begin
          if String.equal c3 "TRF" then nb := 64
        end
        else if String.equal c2 "SY" then begin
          if String.equal c3 "TRF" then nb := (if twostage then 192 else 64)
          else if sname && String.equal c3 "TRD" then nb := 32
          else if sname && String.equal c3 "GST" then nb := 64
        end
        else if cname && String.equal c2 "HE" then begin
          if String.equal c3 "TRF" then nb := (if twostage then 192 else 64)
          else if String.equal c3 "TRD" then nb := 32
          else if String.equal c3 "GST" then nb := 64
        end
        else if unitary then begin
          if gm () && qrset () then nb := 32
        end
        else if String.equal c2 "GB" then begin
          if String.equal c3 "TRF" then nb := (if n4 <= 64 then 1 else 32)
        end
        else if String.equal c2 "PB" then begin
          if String.equal c3 "TRF" then nb := (if n2 <= 64 then 1 else 32)
        end
        else if String.equal c2 "TR" then begin
          if String.equal c3 "TRI" then nb := 64
          else if String.equal c3 "EVC" then nb := 64
          else if String.equal c3 "SYL" then
            nb :=
              (if sname then min (max 48 (min n1 n2 * 16 / 100)) 240
               else min (max 24 (min n1 n2 * 8 / 100)) 80)
        end
        else if String.equal c2 "LA" then begin
          if String.equal c3 "UUM" then nb := 64
          else if String.equal c3 "TRS" then nb := 32
        end
        else if sname && String.equal c2 "ST" then begin
          if String.equal c3 "EBZ" then nb := 1
        end
        else if String.equal c2 "GG" then begin
          nb := 32;
          if String.equal c3 "HD3" then nb := 32
        end;
        !nb
      end
      else if ispec = 2 then begin
        let nbmin = ref 2 in
        if String.equal c2 "GE" then begin
          if factored () then nbmin := 2
          else if String.equal c3 "HRD" then nbmin := 2
          else if String.equal c3 "BRD" then nbmin := 2
          else if String.equal c3 "TRI" then nbmin := 2
          else if qp3rk () then nbmin := 2
        end
        else if String.equal c2 "SY" then begin
          if String.equal c3 "TRF" then nbmin := 8
          else if sname && String.equal c3 "TRD" then nbmin := 2
        end
        else if cname && String.equal c2 "HE" then begin
          if String.equal c3 "TRD" then nbmin := 2
        end
        else if unitary then begin
          if gm () && qrset () then nbmin := 2
        end
        else if String.equal c2 "GG" then begin
          nbmin := 2;
          if String.equal c3 "HD3" then nbmin := 2
        end;
        !nbmin
      end
      else begin
        let nx = ref 0 in
        if String.equal c2 "GE" then begin
          if factored () then nx := 128
          else if String.equal c3 "HRD" then nx := 128
          else if String.equal c3 "BRD" then nx := 128
          else if qp3rk () then nx := 128
        end
        else if String.equal c2 "SY" then begin
          if sname && String.equal c3 "TRD" then nx := 32
        end
        else if cname && String.equal c2 "HE" then begin
          if String.equal c3 "TRD" then nx := 32
        end
        else if unitary then begin
          if Char.equal c3.[0] 'G' && qrset () then nx := 128
        end
        else if String.equal c2 "GG" then begin
          nx := 128;
          if String.equal c3 "HD3" then nx := 128
        end
        else if String.equal c2 "LA" then begin
          if String.equal c3 "RFT" then nx := 64
        end;
        !nx
      end
    end
  end
  else if ispec = 4 then 6
  else if ispec = 5 then 2
  else if ispec = 6 then int_of_float (float_of_int (min n1 n2) *. 1.6)
  else if ispec = 7 then 1
  else if ispec = 8 then 50
  else if ispec = 9 then 25
  else if ispec = 10 then ieeeck 1 0.0 1.0
  else if ispec = 11 then ieeeck 0 0.0 1.0
  else if ispec >= 12 && ispec <= 17 then iparmq ispec name opts n1 n2 n3 n4
  else -1

module Raw = struct
  let xerbla = xerbla
  let dlamch = dlamch
  let dlaisnan = dlaisnan
  let disnan = disnan
  let dlapy2 = dlapy2
  let dlapy3 = dlapy3
  let dlassq = dlassq
  let dlaswp = dlaswp
  let dlaset = dlaset
  let dlacpy = dlacpy
  let dlange = dlange
  let ieeeck = ieeeck
  let iparmq = iparmq
  let ilaenv = ilaenv
end

(* ---- The OCaml surface ---- *)

let lamch = dlamch
let isnan = disnan
let lapy2 = dlapy2
let lapy3 = dlapy3

let swap_rows ?(reverse = false) a ipiv =
  let rows = Mat.rows a in
  let k = Array.length ipiv in
  if k > rows then
    invalid_arg "Laux.swap_rows: more pivots than the matrix has rows";
  Array.iter
    (fun p ->
      if p < 0 || p >= rows then
        invalid_arg "Laux.swap_rows: pivot is not a row of the matrix")
    ipiv;
  dlaswp (Mat.cols a) (Mat.data a) 0 rows 0 (k - 1) ipiv 0
    (if reverse then -1 else 1)

let fill ?(part = Full) ?diag a alpha =
  let beta = match diag with Some b -> b | None -> alpha in
  dlaset part (Mat.rows a) (Mat.cols a) alpha beta (Mat.data a) 0 (Mat.rows a)

let copy ?(part = Full) src dst =
  Mat.same_rows "Laux.copy" src dst;
  if Mat.cols src <> Mat.cols dst then
    invalid_arg "Laux.copy: operands disagree on the number of columns";
  dlacpy part (Mat.rows src) (Mat.cols src) (Mat.data src) 0 (Mat.rows src)
    (Mat.data dst) 0 (Mat.rows dst)

let norm kind a =
  let work =
    match kind with
    | Inf_norm -> Array.make (max 1 (Mat.rows a)) 0.0
    | Max_abs | One_norm | Frobenius -> [||]
  in
  dlange kind (Mat.rows a) (Mat.cols a) (Mat.data a) 0 (Mat.rows a) work
