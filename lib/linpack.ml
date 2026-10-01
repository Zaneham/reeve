(* Reeve, Copyright 2026 Zane Hambly.
   LINPACK, ported from the SLATEC double precision routines. Matrices are
   column major and an element (i, j) zero based lives at i + j * lda. *)

let dgefa a lda n ipvt =
  let info = ref 0 in
  let nm1 = n - 1 in
  if nm1 >= 1 then
    for k = 1 to nm1 do
      let kp1 = k + 1 in
      let l = Blas.idamax (n - k + 1) a (k - 1 + ((k - 1) * lda)) 1 + k in
      ipvt.(k - 1) <- l - 1;
      if a.(l - 1 + ((k - 1) * lda)) = 0.0 then info := k
      else begin
        if l <> k then begin
          let t = a.(l - 1 + ((k - 1) * lda)) in
          a.(l - 1 + ((k - 1) * lda)) <- a.(k - 1 + ((k - 1) * lda));
          a.(k - 1 + ((k - 1) * lda)) <- t
        end;
        let t = -1.0 /. a.(k - 1 + ((k - 1) * lda)) in
        Blas.dscal (n - k) t a (k + ((k - 1) * lda)) 1;
        for j = kp1 to n do
          let t = a.(l - 1 + ((j - 1) * lda)) in
          if l <> k then begin
            a.(l - 1 + ((j - 1) * lda)) <- a.(k - 1 + ((j - 1) * lda));
            a.(k - 1 + ((j - 1) * lda)) <- t
          end;
          Blas.daxpy (n - k) t a (k + ((k - 1) * lda)) 1 a
            (k + ((j - 1) * lda)) 1
        done
      end
    done;
  if n >= 1 then begin
    ipvt.(n - 1) <- n - 1;
    if a.(n - 1 + ((n - 1) * lda)) = 0.0 then info := n
  end;
  !info

let dgesl a lda n ipvt b job =
  let nm1 = n - 1 in
  if job <> 0 then begin
    for k = 1 to n do
      let t = Blas.ddot (k - 1) a ((k - 1) * lda) 1 b 0 1 in
      b.(k - 1) <- (b.(k - 1) -. t) /. a.(k - 1 + ((k - 1) * lda))
    done;
    if nm1 >= 1 then
      for kb = 1 to nm1 do
        let k = n - kb in
        b.(k - 1) <-
          b.(k - 1) +. Blas.ddot (n - k) a (k + ((k - 1) * lda)) 1 b k 1;
        let l = ipvt.(k - 1) + 1 in
        if l <> k then begin
          let t = b.(l - 1) in
          b.(l - 1) <- b.(k - 1);
          b.(k - 1) <- t
        end
      done
  end
  else begin
    if nm1 >= 1 then
      for k = 1 to nm1 do
        let l = ipvt.(k - 1) + 1 in
        let t = b.(l - 1) in
        if l <> k then begin
          b.(l - 1) <- b.(k - 1);
          b.(k - 1) <- t
        end;
        Blas.daxpy (n - k) t a (k + ((k - 1) * lda)) 1 b k 1
      done;
    for kb = 1 to n do
      let k = n + 1 - kb in
      b.(k - 1) <- b.(k - 1) /. a.(k - 1 + ((k - 1) * lda));
      let t = -.b.(k - 1) in
      Blas.daxpy (k - 1) t a ((k - 1) * lda) 1 b 0 1
    done
  end

let dgeco a lda n ipvt z =
  let anorm = ref 0.0 in
  for j = 1 to n do
    anorm := Float.max !anorm (Blas.dasum n a ((j - 1) * lda) 1)
  done;
  ignore (dgefa a lda n ipvt);
  let ek = ref 1.0 in
  for j = 1 to n do
    z.(j - 1) <- 0.0
  done;
  for k = 1 to n do
    if z.(k - 1) <> 0.0 then ek := Float.copy_sign !ek (-.z.(k - 1));
    let akk = a.(k - 1 + ((k - 1) * lda)) in
    if Float.abs (!ek -. z.(k - 1)) > Float.abs akk then begin
      let s = Float.abs akk /. Float.abs (!ek -. z.(k - 1)) in
      Blas.dscal n s z 0 1;
      ek := s *. !ek
    end;
    let wk = ref (!ek -. z.(k - 1)) and wkm = ref (-. !ek -. z.(k - 1)) in
    let s = ref (Float.abs !wk) and sm = ref (Float.abs !wkm) in
    if akk = 0.0 then begin
      wk := 1.0;
      wkm := 1.0
    end
    else begin
      wk := !wk /. akk;
      wkm := !wkm /. akk
    end;
    let kp1 = k + 1 in
    if kp1 <= n then begin
      for j = kp1 to n do
        let akj = a.(k - 1 + ((j - 1) * lda)) in
        sm := !sm +. Float.abs (z.(j - 1) +. (!wkm *. akj));
        z.(j - 1) <- z.(j - 1) +. (!wk *. akj);
        s := !s +. Float.abs z.(j - 1)
      done;
      if !s < !sm then begin
        let t = !wkm -. !wk in
        wk := !wkm;
        for j = kp1 to n do
          z.(j - 1) <- z.(j - 1) +. (t *. a.(k - 1 + ((j - 1) * lda)))
        done
      end
    end;
    z.(k - 1) <- !wk
  done;
  Blas.dscal n (1.0 /. Blas.dasum n z 0 1) z 0 1;
  for kb = 1 to n do
    let k = n + 1 - kb in
    if k < n then
      z.(k - 1) <-
        z.(k - 1) +. Blas.ddot (n - k) a (k + ((k - 1) * lda)) 1 z k 1;
    if Float.abs z.(k - 1) > 1.0 then
      Blas.dscal n (1.0 /. Float.abs z.(k - 1)) z 0 1;
    let l = ipvt.(k - 1) + 1 in
    let t = z.(l - 1) in
    z.(l - 1) <- z.(k - 1);
    z.(k - 1) <- t
  done;
  Blas.dscal n (1.0 /. Blas.dasum n z 0 1) z 0 1;
  let ynorm = ref 1.0 in
  for k = 1 to n do
    let l = ipvt.(k - 1) + 1 in
    let t = z.(l - 1) in
    z.(l - 1) <- z.(k - 1);
    z.(k - 1) <- t;
    if k < n then Blas.daxpy (n - k) t a (k + ((k - 1) * lda)) 1 z k 1;
    if Float.abs z.(k - 1) > 1.0 then begin
      let s = 1.0 /. Float.abs z.(k - 1) in
      Blas.dscal n s z 0 1;
      ynorm := s *. !ynorm
    end
  done;
  let s = 1.0 /. Blas.dasum n z 0 1 in
  Blas.dscal n s z 0 1;
  ynorm := s *. !ynorm;
  for kb = 1 to n do
    let k = n + 1 - kb in
    let akk = a.(k - 1 + ((k - 1) * lda)) in
    if Float.abs z.(k - 1) > Float.abs akk then begin
      let s = Float.abs akk /. Float.abs z.(k - 1) in
      Blas.dscal n s z 0 1;
      ynorm := s *. !ynorm
    end;
    if akk <> 0.0 then z.(k - 1) <- z.(k - 1) /. akk;
    if akk = 0.0 then z.(k - 1) <- 1.0;
    let t = -.z.(k - 1) in
    Blas.daxpy (k - 1) t a ((k - 1) * lda) 1 z 0 1
  done;
  let s = 1.0 /. Blas.dasum n z 0 1 in
  Blas.dscal n s z 0 1;
  ynorm := s *. !ynorm;
  if !anorm <> 0.0 then !ynorm /. !anorm else 0.0

let dgbfa abd lda n ml mu ipvt =
  let m = ml + mu + 1 in
  let info = ref 0 in
  let j0 = mu + 2 in
  let j1 = min n m - 1 in
  if j1 >= j0 then
    for jz = j0 to j1 do
      let i0 = m + 1 - jz in
      for i = i0 to ml do
        abd.(i - 1 + ((jz - 1) * lda)) <- 0.0
      done
    done;
  let jz = ref j1 in
  let ju = ref 0 in
  let nm1 = n - 1 in
  if nm1 >= 1 then
    for k = 1 to nm1 do
      let kp1 = k + 1 in
      jz := !jz + 1;
      if !jz <= n && ml >= 1 then
        for i = 1 to ml do
          abd.(i - 1 + ((!jz - 1) * lda)) <- 0.0
        done;
      let lm = min ml (n - k) in
      let l = ref (Blas.idamax (lm + 1) abd (m - 1 + ((k - 1) * lda)) 1 + m) in
      ipvt.(k - 1) <- !l + k - m - 1;
      if abd.(!l - 1 + ((k - 1) * lda)) = 0.0 then info := k
      else begin
        if !l <> m then begin
          let t = abd.(!l - 1 + ((k - 1) * lda)) in
          abd.(!l - 1 + ((k - 1) * lda)) <- abd.(m - 1 + ((k - 1) * lda));
          abd.(m - 1 + ((k - 1) * lda)) <- t
        end;
        let t = -1.0 /. abd.(m - 1 + ((k - 1) * lda)) in
        Blas.dscal lm t abd (m + ((k - 1) * lda)) 1;
        ju := min (max !ju (mu + !l + k - m)) n;
        let mm = ref m in
        if !ju >= kp1 then
          for j = kp1 to !ju do
            l := !l - 1;
            mm := !mm - 1;
            let t = abd.(!l - 1 + ((j - 1) * lda)) in
            if !l <> !mm then begin
              abd.(!l - 1 + ((j - 1) * lda)) <- abd.(!mm - 1 + ((j - 1) * lda));
              abd.(!mm - 1 + ((j - 1) * lda)) <- t
            end;
            Blas.daxpy lm t abd (m + ((k - 1) * lda)) 1 abd
              (!mm + ((j - 1) * lda)) 1
          done
      end
    done;
  if n >= 1 then begin
    ipvt.(n - 1) <- n - 1;
    if abd.(m - 1 + ((n - 1) * lda)) = 0.0 then info := n
  end;
  !info

let dgbsl abd lda n ml mu ipvt b job =
  let m = mu + ml + 1 in
  let nm1 = n - 1 in
  if job <> 0 then begin
    for k = 1 to n do
      let lm = min k m - 1 in
      let la = m - lm in
      let lb = k - lm in
      let t = Blas.ddot lm abd (la - 1 + ((k - 1) * lda)) 1 b (lb - 1) 1 in
      b.(k - 1) <- (b.(k - 1) -. t) /. abd.(m - 1 + ((k - 1) * lda))
    done;
    if ml <> 0 && nm1 >= 1 then
      for kb = 1 to nm1 do
        let k = n - kb in
        let lm = min ml (n - k) in
        b.(k - 1) <-
          b.(k - 1) +. Blas.ddot lm abd (m + ((k - 1) * lda)) 1 b k 1;
        let l = ipvt.(k - 1) + 1 in
        if l <> k then begin
          let t = b.(l - 1) in
          b.(l - 1) <- b.(k - 1);
          b.(k - 1) <- t
        end
      done
  end
  else begin
    if ml <> 0 && nm1 >= 1 then
      for k = 1 to nm1 do
        let lm = min ml (n - k) in
        let l = ipvt.(k - 1) + 1 in
        let t = b.(l - 1) in
        if l <> k then begin
          b.(l - 1) <- b.(k - 1);
          b.(k - 1) <- t
        end;
        Blas.daxpy lm t abd (m + ((k - 1) * lda)) 1 b k 1
      done;
    for kb = 1 to n do
      let k = n + 1 - kb in
      b.(k - 1) <- b.(k - 1) /. abd.(m - 1 + ((k - 1) * lda));
      let lm = min k m - 1 in
      let la = m - lm in
      let lb = k - lm in
      let t = -.b.(k - 1) in
      Blas.daxpy lm t abd (la - 1 + ((k - 1) * lda)) 1 b (lb - 1) 1
    done
  end

let dpofa a lda n =
  let info = ref 0 in
  let j = ref 1 in
  let bad = ref false in
  while (not !bad) && !j <= n do
    info := !j;
    let s = ref 0.0 in
    let jm1 = !j - 1 in
    if jm1 >= 1 then
      for k = 1 to jm1 do
        let t =
          (a.(k - 1 + ((!j - 1) * lda))
          -. Blas.ddot (k - 1) a ((k - 1) * lda) 1 a ((!j - 1) * lda) 1)
          /. a.(k - 1 + ((k - 1) * lda))
        in
        a.(k - 1 + ((!j - 1) * lda)) <- t;
        s := !s +. (t *. t)
      done;
    let d = a.(!j - 1 + ((!j - 1) * lda)) -. !s in
    if d <= 0.0 then bad := true
    else begin
      a.(!j - 1 + ((!j - 1) * lda)) <- sqrt d;
      incr j
    end
  done;
  if !bad then !info else 0

let dposl a lda n b =
  for k = 1 to n do
    let t = Blas.ddot (k - 1) a ((k - 1) * lda) 1 b 0 1 in
    b.(k - 1) <- (b.(k - 1) -. t) /. a.(k - 1 + ((k - 1) * lda))
  done;
  for kb = 1 to n do
    let k = n + 1 - kb in
    b.(k - 1) <- b.(k - 1) /. a.(k - 1 + ((k - 1) * lda));
    let t = -.b.(k - 1) in
    Blas.daxpy (k - 1) t a ((k - 1) * lda) 1 b 0 1
  done

(* These two aren't in the library at all, only as dpbfa_local and dpbsl_local
   inside the extended tests, and that dpbsl gathers where the real one
   scatters. *)
let dpbfa abd lda n m =
  let info = ref 0 in
  let j = ref 1 in
  let bad = ref false in
  while (not !bad) && !j <= n do
    let s = ref 0.0 in
    let ik = ref (m + 1) in
    let jk = ref (max (!j - m) 1) in
    let mu = max (m + 2 - !j) 1 in
    if m >= mu then
      for k = mu to m do
        let t =
          (abd.(k - 1 + ((!j - 1) * lda))
          -. Blas.ddot (k - mu) abd (!ik - 1 + ((!jk - 1) * lda)) 1 abd
               (mu - 1 + ((!j - 1) * lda)) 1)
          /. abd.(m + ((!jk - 1) * lda))
        in
        abd.(k - 1 + ((!j - 1) * lda)) <- t;
        s := !s +. (t *. t);
        ik := !ik - 1;
        jk := !jk + 1
      done;
    let d = abd.(m + ((!j - 1) * lda)) -. !s in
    if d <= 0.0 then begin
      info := !j;
      bad := true
    end
    else begin
      abd.(m + ((!j - 1) * lda)) <- sqrt d;
      incr j
    end
  done;
  if !bad then !info else 0

let dpbsl abd lda n m b =
  for k = 1 to n do
    let lm = min (k - 1) m in
    let la = m + 1 - lm in
    let lb = k - lm in
    let t = Blas.ddot lm abd (la - 1 + ((k - 1) * lda)) 1 b (lb - 1) 1 in
    b.(k - 1) <- (b.(k - 1) -. t) /. abd.(m + ((k - 1) * lda))
  done;
  for kb = 1 to n do
    let k = n + 1 - kb in
    let lm = min m (n - k) in
    let t = ref 0.0 in
    for l = 1 to lm do
      t := !t +. (abd.(m - l + ((k + l - 1) * lda)) *. b.(k + l - 1))
    done;
    b.(k - 1) <- (b.(k - 1) -. !t) /. abd.(m + ((k - 1) * lda))
  done
