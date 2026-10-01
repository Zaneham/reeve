(* Reeve, Copyright 2026 Zane Hambly.
   LAPACK's dense linear system routines, the LU and Cholesky paths, ported
   from the reference Fortran. Matrices carry a base offset and a leading
   dimension as in {!Blasmat}, and pivot vectors hold zero based row
   indices. *)

open Blasmat

let arg name pos = invalid_arg (name ^ ": argument " ^ string_of_int pos)
let is_upper u = match u with Upper -> true | Lower -> false
let sfmin = Float.min_float
let nb_getrf = 64
let nb_potrf = 64

let dgetf2 m n a ao lda ipiv ipo =
  if m < 0 then arg "dgetf2" 1
  else if n < 0 then arg "dgetf2" 2
  else if lda < max 1 m then arg "dgetf2" 4
  else if m = 0 || n = 0 then 0
  else begin
    let info = ref 0 in
    let mn = min m n in
    for j = 1 to mn do
      let d = ao + j - 1 + ((j - 1) * lda) in
      let jp = j - 1 + Blas.idamax (m - j + 1) a d 1 in
      ipiv.(ipo + j - 1) <- jp;
      if a.(ao + jp + ((j - 1) * lda)) <> 0.0 then begin
        if jp <> j - 1 then
          Blas.dswap n a (ao + j - 1) lda a (ao + jp) lda;
        if j < m then
          if Float.abs a.(d) >= sfmin then
            Blas.dscal (m - j) (1.0 /. a.(d)) a (d + 1) 1
          else
            for i = 1 to m - j do
              a.(d + i) <- a.(d + i) /. a.(d)
            done
      end
      else if !info = 0 then info := j;
      if j < mn then
        dger (m - j) (n - j) (-1.0) a (d + 1) 1 a (d + lda) lda a
          (d + 1 + lda) lda
    done;
    !info
  end

let rec dgetrf2 m n a ao lda ipiv ipo =
  if m < 0 then arg "dgetrf2" 1
  else if n < 0 then arg "dgetrf2" 2
  else if lda < max 1 m then arg "dgetrf2" 4
  else if m = 0 || n = 0 then 0
  else if m = 1 then begin
    ipiv.(ipo) <- 0;
    if a.(ao) = 0.0 then 1 else 0
  end
  else if n = 1 then begin
    let i = Blas.idamax m a ao 1 in
    ipiv.(ipo) <- i;
    if a.(ao + i) <> 0.0 then begin
      if i <> 0 then begin
        let t = a.(ao) in
        a.(ao) <- a.(ao + i);
        a.(ao + i) <- t
      end;
      if Float.abs a.(ao) >= sfmin then
        Blas.dscal (m - 1) (1.0 /. a.(ao)) a (ao + 1) 1
      else
        for k = 1 to m - 1 do
          a.(ao + k) <- a.(ao + k) /. a.(ao)
        done;
      0
    end
    else 1
  end
  else begin
    let n1 = min m n / 2 in
    let n2 = n - n1 in
    let info = ref 0 in
    let left = dgetrf2 m n1 a ao lda ipiv ipo in
    if !info = 0 && left > 0 then info := left;
    Laux.dlaswp n2 a (ao + (n1 * lda)) lda 0 (n1 - 1) ipiv ipo 1;
    dtrsm Left Lower No_trans Unit n1 n2 1.0 a ao lda a
      (ao + (n1 * lda)) lda;
    dgemm No_trans No_trans (m - n1) n2 n1 (-1.0) a (ao + n1) lda a
      (ao + (n1 * lda)) lda 1.0 a (ao + n1 + (n1 * lda)) lda;
    let right =
      dgetrf2 (m - n1) n2 a (ao + n1 + (n1 * lda)) lda ipiv (ipo + n1)
    in
    if !info = 0 && right > 0 then info := right + n1;
    for i = n1 + 1 to min m n do
      ipiv.(ipo + i - 1) <- ipiv.(ipo + i - 1) + n1
    done;
    Laux.dlaswp n1 a ao lda n1 (min m n - 1) ipiv ipo 1;
    !info
  end

let dgetrf m n a ao lda ipiv ipo =
  if m < 0 then arg "dgetrf" 1
  else if n < 0 then arg "dgetrf" 2
  else if lda < max 1 m then arg "dgetrf" 4
  else if m = 0 || n = 0 then 0
  else begin
    let mn = min m n in
    if nb_getrf <= 1 || nb_getrf >= mn then dgetrf2 m n a ao lda ipiv ipo
    else begin
      let info = ref 0 in
      let j = ref 1 in
      while !j <= mn do
        let jb = min (mn - !j + 1) nb_getrf in
        let d = ao + !j - 1 + ((!j - 1) * lda) in
        let iinfo =
          dgetrf2 (m - !j + 1) jb a d lda ipiv (ipo + !j - 1)
        in
        if !info = 0 && iinfo > 0 then info := iinfo + !j - 1;
        for i = !j to min m (!j + jb - 1) do
          ipiv.(ipo + i - 1) <- !j - 1 + ipiv.(ipo + i - 1)
        done;
        Laux.dlaswp (!j - 1) a ao lda (!j - 1) (!j + jb - 2) ipiv ipo 1;
        if !j + jb <= n then begin
          let r = ao + ((!j + jb - 1) * lda) in
          Laux.dlaswp (n - !j - jb + 1) a r lda (!j - 1) (!j + jb - 2)
            ipiv ipo 1;
          dtrsm Left Lower No_trans Unit jb (n - !j - jb + 1) 1.0 a d lda a
            (r + !j - 1) lda;
          if !j + jb <= m then
            dgemm No_trans No_trans (m - !j - jb + 1) (n - !j - jb + 1) jb
              (-1.0) a (d + jb) lda a (r + !j - 1) lda 1.0 a
              (r + !j + jb - 1) lda
        end;
        j := !j + nb_getrf
      done;
      !info
    end
  end

let dgetrs trans n nrhs a ao lda ipiv ipo b bo ldb =
  if n < 0 then arg "dgetrs" 2
  else if nrhs < 0 then arg "dgetrs" 3
  else if lda < max 1 n then arg "dgetrs" 5
  else if ldb < max 1 n then arg "dgetrs" 8
  else if n = 0 || nrhs = 0 then 0
  else begin
    (match trans with
    | No_trans ->
      Laux.dlaswp nrhs b bo ldb 0 (n - 1) ipiv ipo 1;
      dtrsm Left Lower No_trans Unit n nrhs 1.0 a ao lda b bo ldb;
      dtrsm Left Upper No_trans Non_unit n nrhs 1.0 a ao lda b bo ldb
    | Trans | Conj_trans ->
      dtrsm Left Upper Trans Non_unit n nrhs 1.0 a ao lda b bo ldb;
      dtrsm Left Lower Trans Unit n nrhs 1.0 a ao lda b bo ldb;
      Laux.dlaswp nrhs b bo ldb 0 (n - 1) ipiv ipo (-1));
    0
  end

let dgesv n nrhs a ao lda ipiv ipo b bo ldb =
  if n < 0 then arg "dgesv" 1
  else if nrhs < 0 then arg "dgesv" 2
  else if lda < max 1 n then arg "dgesv" 4
  else if ldb < max 1 n then arg "dgesv" 7
  else begin
    let info = dgetrf n n a ao lda ipiv ipo in
    if info = 0 then dgetrs No_trans n nrhs a ao lda ipiv ipo b bo ldb
    else info
  end

let dpotf2 uplo n a ao lda =
  if n < 0 then arg "dpotf2" 2
  else if lda < max 1 n then arg "dpotf2" 4
  else if n = 0 then 0
  else begin
    let info = ref 0 in
    let j = ref 1 in
    while !info = 0 && !j <= n do
      let d = ao + !j - 1 + ((!j - 1) * lda) in
      let col = ao + ((!j - 1) * lda) in
      let row = ao + !j - 1 in
      let ajj =
        if is_upper uplo then
          a.(d) -. Blas.ddot (!j - 1) a col 1 a col 1
        else a.(d) -. Blas.ddot (!j - 1) a row lda a row lda
      in
      if ajj <= 0.0 || Float.is_nan ajj then begin
        a.(d) <- ajj;
        info := !j
      end
      else begin
        let r = sqrt ajj in
        a.(d) <- r;
        if !j < n then
          if is_upper uplo then begin
            dgemv Trans (!j - 1) (n - !j) (-1.0) a (ao + (!j * lda)) lda a col
              1 1.0 a (d + lda) lda;
            Blas.dscal (n - !j) (1.0 /. r) a (d + lda) lda
          end
          else begin
            dgemv No_trans (n - !j) (!j - 1) (-1.0) a (row + 1) lda a row lda
              1.0 a (d + 1) 1;
            Blas.dscal (n - !j) (1.0 /. r) a (d + 1) 1
          end;
        incr j
      end
    done;
    !info
  end

let rec dpotrf2 uplo n a ao lda =
  if n < 0 then arg "dpotrf2" 2
  else if lda < max 1 n then arg "dpotrf2" 4
  else if n = 0 then 0
  else if n = 1 then
    if a.(ao) <= 0.0 || Float.is_nan a.(ao) then 1
    else begin
      a.(ao) <- sqrt a.(ao);
      0
    end
  else begin
    let n1 = n / 2 in
    let n2 = n - n1 in
    let first = dpotrf2 uplo n1 a ao lda in
    if first <> 0 then first
    else begin
      let d = ao + n1 + (n1 * lda) in
      if is_upper uplo then begin
        dtrsm Left Upper Trans Non_unit n1 n2 1.0 a ao lda a
          (ao + (n1 * lda)) lda;
        dsyrk uplo Trans n2 n1 (-1.0) a (ao + (n1 * lda)) lda 1.0 a d lda
      end
      else begin
        dtrsm Right Lower Trans Non_unit n2 n1 1.0 a ao lda a (ao + n1) lda;
        dsyrk uplo No_trans n2 n1 (-1.0) a (ao + n1) lda 1.0 a d lda
      end;
      let second = dpotrf2 uplo n2 a d lda in
      if second <> 0 then second + n1 else 0
    end
  end

let dpotrf uplo n a ao lda =
  if n < 0 then arg "dpotrf" 2
  else if lda < max 1 n then arg "dpotrf" 4
  else if n = 0 then 0
  else if nb_potrf <= 1 || nb_potrf >= n then dpotrf2 uplo n a ao lda
  else begin
    let info = ref 0 in
    let j = ref 1 in
    while !info = 0 && !j <= n do
      let jb = min nb_potrf (n - !j + 1) in
      let d = ao + !j - 1 + ((!j - 1) * lda) in
      if is_upper uplo then
        dsyrk Upper Trans jb (!j - 1) (-1.0) a (ao + ((!j - 1) * lda)) lda 1.0
          a d lda
      else
        dsyrk Lower No_trans jb (!j - 1) (-1.0) a (ao + !j - 1) lda 1.0 a d
          lda;
      let iinfo = dpotrf2 uplo jb a d lda in
      if iinfo <> 0 then info := iinfo + !j - 1
      else begin
        if !j + jb <= n then
          if is_upper uplo then begin
            let r = ao + !j - 1 + ((!j + jb - 1) * lda) in
            dgemm Trans No_trans jb (n - !j - jb + 1) (!j - 1) (-1.0) a
              (ao + ((!j - 1) * lda)) lda a (ao + ((!j + jb - 1) * lda)) lda
              1.0 a r lda;
            dtrsm Left Upper Trans Non_unit jb (n - !j - jb + 1) 1.0 a d lda a
              r lda
          end
          else begin
            let r = ao + !j + jb - 1 + ((!j - 1) * lda) in
            dgemm No_trans Trans (n - !j - jb + 1) jb (!j - 1) (-1.0) a
              (ao + !j + jb - 1) lda a (ao + !j - 1) lda 1.0 a r lda;
            dtrsm Right Lower Trans Non_unit (n - !j - jb + 1) jb 1.0 a d lda
              a r lda
          end;
        j := !j + nb_potrf
      end
    done;
    !info
  end

let dpotrs uplo n nrhs a ao lda b bo ldb =
  if n < 0 then arg "dpotrs" 2
  else if nrhs < 0 then arg "dpotrs" 3
  else if lda < max 1 n then arg "dpotrs" 5
  else if ldb < max 1 n then arg "dpotrs" 7
  else if n = 0 || nrhs = 0 then 0
  else begin
    if is_upper uplo then begin
      dtrsm Left Upper Trans Non_unit n nrhs 1.0 a ao lda b bo ldb;
      dtrsm Left Upper No_trans Non_unit n nrhs 1.0 a ao lda b bo ldb
    end
    else begin
      dtrsm Left Lower No_trans Non_unit n nrhs 1.0 a ao lda b bo ldb;
      dtrsm Left Lower Trans Non_unit n nrhs 1.0 a ao lda b bo ldb
    end;
    0
  end

let dposv uplo n nrhs a ao lda b bo ldb =
  if n < 0 then arg "dposv" 2
  else if nrhs < 0 then arg "dposv" 3
  else if lda < max 1 n then arg "dposv" 5
  else if ldb < max 1 n then arg "dposv" 7
  else begin
    let info = dpotrf uplo n a ao lda in
    if info = 0 then dpotrs uplo n nrhs a ao lda b bo ldb else info
  end
