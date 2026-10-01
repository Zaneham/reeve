(* Reeve, Copyright 2026 Zane Hambly.
   The double precision data handling routines, ported from the SLATEC
   Fortran. Permutation vectors hold zero based indices, as pivot vectors do
   elsewhere in the library, and the Fortran's signed KFLAG is a variant. *)

type kflag = Increasing | Increasing_carry | Decreasing | Decreasing_carry

let copy n x xo y yo =
  for p = 0 to n - 1 do
    y.(yo + p) <- x.(xo + p)
  done

let carries kflag =
  match kflag with
  | Increasing | Decreasing -> false
  | Increasing_carry | Decreasing_carry -> true

let descends kflag =
  match kflag with
  | Increasing | Increasing_carry -> false
  | Decreasing | Decreasing_carry -> true

(* ---- dsort ---- *)

(* Label 90 is unreachable. The guard at dsort.f:179 kicks every segment of two
   or more back to the partitioner, so the insertion sort only ever gets handed
   a single element. Singleton's cutoff was 11. *)
let dsort dx dy n kflag =
  let nn = n in
  if nn < 1 then Laux.xerbla "dsort" 3;
  let carry = carries kflag in
  let down = descends kflag in
  if down then
    for p = 1 to nn do
      dx.(p - 1) <- -.dx.(p - 1)
    done;
  let il = Array.make 21 0 in
  let iu = Array.make 21 0 in
  let m = ref 1 in
  let i = ref 1 in
  let j = ref nn in
  let r = ref 0.375 in
  let ij = ref 1 in
  let k = ref 1 in
  let l = ref 1 in
  let t = ref 0.0 in
  let ty = ref 0.0 in
  let rec main () =
    if !i = !j then pop ()
    else begin
      if !r <= 0.5898437 then r := !r +. 3.90625e-2 else r := !r -. 0.21875;
      pivot ()
    end
  and pivot () =
    k := !i;
    ij := !i + int_of_float (float_of_int (!j - !i) *. !r);
    t := dx.(!ij - 1);
    if carry then ty := dy.(!ij - 1);
    if dx.(!i - 1) > !t then begin
      dx.(!ij - 1) <- dx.(!i - 1);
      dx.(!i - 1) <- !t;
      t := dx.(!ij - 1);
      if carry then begin
        dy.(!ij - 1) <- dy.(!i - 1);
        dy.(!i - 1) <- !ty;
        ty := dy.(!ij - 1)
      end
    end;
    l := !j;
    if dx.(!j - 1) < !t then begin
      dx.(!ij - 1) <- dx.(!j - 1);
      dx.(!j - 1) <- !t;
      t := dx.(!ij - 1);
      if carry then begin
        dy.(!ij - 1) <- dy.(!j - 1);
        dy.(!j - 1) <- !ty;
        ty := dy.(!ij - 1)
      end;
      if dx.(!i - 1) > !t then begin
        dx.(!ij - 1) <- dx.(!i - 1);
        dx.(!i - 1) <- !t;
        t := dx.(!ij - 1);
        if carry then begin
          dy.(!ij - 1) <- dy.(!i - 1);
          dy.(!i - 1) <- !ty;
          ty := dy.(!ij - 1)
        end
      end
    end;
    lower ()
  and lower () =
    l := !l - 1;
    if dx.(!l - 1) > !t then lower () else upper ()
  and upper () =
    k := !k + 1;
    if dx.(!k - 1) < !t then upper ()
    else if !k <= !l then begin
      let tt = dx.(!l - 1) in
      dx.(!l - 1) <- dx.(!k - 1);
      dx.(!k - 1) <- tt;
      if carry then begin
        let tty = dy.(!l - 1) in
        dy.(!l - 1) <- dy.(!k - 1);
        dy.(!k - 1) <- tty
      end;
      lower ()
    end
    else begin
      if !l - !i > !j - !k then begin
        il.(!m - 1) <- !i;
        iu.(!m - 1) <- !l;
        i := !k;
        m := !m + 1
      end
      else begin
        il.(!m - 1) <- !k;
        iu.(!m - 1) <- !j;
        j := !l;
        m := !m + 1
      end;
      decide ()
    end
  and pop () =
    m := !m - 1;
    if !m <> 0 then begin
      i := il.(!m - 1);
      j := iu.(!m - 1);
      decide ()
    end
  and decide () =
    if !j - !i >= 1 then pivot ()
    else if !i = 1 then main ()
    else begin
      i := !i - 1;
      insert ()
    end
  and insert () =
    i := !i + 1;
    if !i = !j then pop ()
    else begin
      t := dx.(!i);
      if carry then ty := dy.(!i);
      if dx.(!i - 1) <= !t then insert ()
      else begin
        k := !i;
        shift ()
      end
    end
  and shift () =
    dx.(!k) <- dx.(!k - 1);
    if carry then dy.(!k) <- dy.(!k - 1);
    k := !k - 1;
    if !t < dx.(!k - 1) then shift ()
    else begin
      dx.(!k) <- !t;
      if carry then dy.(!k) <- !ty;
      insert ()
    end
  in
  main ();
  if down then
    for p = 1 to nn do
      dx.(p - 1) <- -.dx.(p - 1)
    done

(* ---- dpsort ---- *)

let dpsort dx n iperm kflag =
  let nn = n in
  if nn < 1 then Laux.xerbla "dpsort" 2;
  for p = 1 to nn do
    iperm.(p - 1) <- p
  done;
  if nn > 1 then begin
    let rearrange = carries kflag in
    let down = descends kflag in
    if down then
      for p = 1 to nn do
        dx.(p - 1) <- -.dx.(p - 1)
      done;
    let il = Array.make 21 0 in
    let iu = Array.make 21 0 in
    let m = ref 1 in
    let i = ref 1 in
    let j = ref nn in
    let r = ref 0.375 in
    let ij = ref 1 in
    let k = ref 1 in
    let l = ref 1 in
    let lm = ref 1 in
    let rec main () =
      if !i = !j then pop ()
      else begin
        if !r <= 0.5898437 then r := !r +. 3.90625e-2 else r := !r -. 0.21875;
        pivot ()
      end
    and pivot () =
      k := !i;
      ij := !i + int_of_float (float_of_int (!j - !i) *. !r);
      lm := iperm.(!ij - 1);
      if dx.(iperm.(!i - 1) - 1) > dx.(!lm - 1) then begin
        iperm.(!ij - 1) <- iperm.(!i - 1);
        iperm.(!i - 1) <- !lm;
        lm := iperm.(!ij - 1)
      end;
      l := !j;
      if dx.(iperm.(!j - 1) - 1) < dx.(!lm - 1) then begin
        iperm.(!ij - 1) <- iperm.(!j - 1);
        iperm.(!j - 1) <- !lm;
        lm := iperm.(!ij - 1);
        if dx.(iperm.(!i - 1) - 1) > dx.(!lm - 1) then begin
          iperm.(!ij - 1) <- iperm.(!i - 1);
          iperm.(!i - 1) <- !lm;
          lm := iperm.(!ij - 1)
        end
      end;
      lower ()
    and lower () =
      l := !l - 1;
      if dx.(iperm.(!l - 1) - 1) > dx.(!lm - 1) then lower () else upper ()
    and upper () =
      k := !k + 1;
      if dx.(iperm.(!k - 1) - 1) < dx.(!lm - 1) then upper ()
      else if !k <= !l then begin
        let lmt = iperm.(!l - 1) in
        iperm.(!l - 1) <- iperm.(!k - 1);
        iperm.(!k - 1) <- lmt;
        lower ()
      end
      else begin
        if !l - !i > !j - !k then begin
          il.(!m - 1) <- !i;
          iu.(!m - 1) <- !l;
          i := !k;
          m := !m + 1
        end
        else begin
          il.(!m - 1) <- !k;
          iu.(!m - 1) <- !j;
          j := !l;
          m := !m + 1
        end;
        decide ()
      end
    and pop () =
      m := !m - 1;
      if !m <> 0 then begin
        i := il.(!m - 1);
        j := iu.(!m - 1);
        decide ()
      end
    and decide () =
      if !j - !i >= 1 then pivot ()
      else if !i = 1 then main ()
      else begin
        i := !i - 1;
        insert ()
      end
    and insert () =
      i := !i + 1;
      if !i = !j then pop ()
      else begin
        lm := iperm.(!i);
        if dx.(iperm.(!i - 1) - 1) <= dx.(!lm - 1) then insert ()
        else begin
          k := !i;
          shift ()
        end
      end
    and shift () =
      iperm.(!k) <- iperm.(!k - 1);
      k := !k - 1;
      if dx.(!lm - 1) < dx.(iperm.(!k - 1) - 1) then shift ()
      else begin
        iperm.(!k) <- !lm;
        insert ()
      end
    in
    main ();
    if down then
      for p = 1 to nn do
        dx.(p - 1) <- -.dx.(p - 1)
      done;
    if rearrange then begin
      for istrt = 1 to nn do
        if iperm.(istrt - 1) >= 0 then begin
          let indx = ref istrt in
          let indx0 = ref istrt in
          let temp = dx.(istrt - 1) in
          while iperm.(!indx - 1) > 0 do
            dx.(!indx - 1) <- dx.(iperm.(!indx - 1) - 1);
            indx0 := !indx;
            iperm.(!indx - 1) <- -iperm.(!indx - 1);
            indx := abs iperm.(!indx - 1)
          done;
          dx.(!indx0 - 1) <- temp
        end
      done;
      for p = 1 to nn do
        iperm.(p - 1) <- -iperm.(p - 1)
      done
    end
  end;
  for p = 1 to nn do
    iperm.(p - 1) <- iperm.(p - 1) - 1
  done

(* ---- dpperm ---- *)

let dpperm dx n iperm =
  if n < 1 then Laux.xerbla "dpperm" 2;
  let ip = Array.init n (fun p -> iperm.(p) + 1) in
  for i = 1 to n do
    let indx = abs ip.(i - 1) in
    if indx >= 1 && indx <= n && ip.(indx - 1) > 0 then
      ip.(indx - 1) <- -ip.(indx - 1)
    else Laux.xerbla "dpperm" 3
  done;
  for istrt = 1 to n do
    if ip.(istrt - 1) <= 0 then begin
      let indx = ref istrt in
      let indx0 = ref istrt in
      let dtemp = dx.(istrt - 1) in
      while ip.(!indx - 1) < 0 do
        dx.(!indx - 1) <- dx.(-ip.(!indx - 1) - 1);
        indx0 := !indx;
        ip.(!indx - 1) <- -ip.(!indx - 1);
        indx := ip.(!indx - 1)
      done;
      dx.(!indx0 - 1) <- dtemp
    end
  done

(* ---- d1merg ---- *)

let d1merg tcos i1 m1 i2 m2 i3 =
  if m1 = 0 && m2 = 0 then ()
  else if m1 = 0 then copy m2 tcos i2 tcos i3
  else if m2 = 0 then copy m1 tcos i1 tcos i3
  else begin
    let rec step j1 j2 j3 =
      if tcos.(i1 + j1) <= tcos.(i2 + j2) then begin
        tcos.(i3 + j3) <- tcos.(i1 + j1);
        let j1 = j1 + 1 in
        if j1 > m1 - 1 then copy (m2 - j2) tcos (i2 + j2) tcos (i3 + j3 + 1)
        else step j1 j2 (j3 + 1)
      end
      else begin
        tcos.(i3 + j3) <- tcos.(i2 + j2);
        let j2 = j2 + 1 in
        if j2 > m2 - 1 then copy (m1 - j1) tcos (i1 + j1) tcos (i3 + j3 + 1)
        else step j1 j2 (j3 + 1)
      end
    in
    step 0 0 0
  end
