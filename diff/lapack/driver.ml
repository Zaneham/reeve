open Blas
open Blasmat
open Lapack

let sd = ref 0
let rseed v = sd := v

let rnd () =
  sd := ((!sd * 1103515245) + 12345) mod 2147483648;
  let v = float_of_int ((!sd mod 20001) - 10000) /. 1000.0 in
  if v = 0.0 then 0.25 else v

let buf = Buffer.create (1 lsl 20)

let pr x =
  Buffer.add_string buf (Printf.sprintf "%016LX\n" (Int64.bits_of_float x))

let pin i = Buffer.add_string buf (Printf.sprintf "%d\n" i)

let lbl c i1 i2 i3 i4 i5 i6 =
  Buffer.add_string buf (Printf.sprintf "%d %d %d %d %d %d %d\n" c i1 i2 i3 i4 i5 i6);
  if Buffer.length buf > 1 lsl 22 then begin
    print_string (Buffer.contents buf);
    Buffer.clear buf
  end

let fillb x n = for i = 0 to n - 1 do x.(i) <- rnd () done
let dumpb x n = for i = 0 to n - 1 do pr x.(i) done

let fillspd x nbt off ld n =
  fillb x nbt;
  for j = 1 to n do
    for i = j + 1 to n do
      let v = rnd () in
      x.(off + (i - 1) + ((j - 1) * ld)) <- v;
      x.(off + (j - 1) + ((i - 1) * ld)) <- v
    done
  done;
  for j = 1 to n do
    x.(off + (j - 1) + ((j - 1) * ld)) <-
      (10.0 *. float_of_int n) +. float_of_int j
  done

let tr i = if i = 1 then No_trans else if i = 2 then Trans else Conj_trans
let up i = if i = 1 then Upper else Lower
let sd_ i = if i = 1 then Left else Right
let dg i = if i = 1 then Unit else Non_unit

let dm3 = [| 0; 1; 3 |]
let dm4 = [| 0; 1; 3; 5 |]
let al3 = [| 0.0; 1.0; 2.5 |]
let be3 = [| 0.0; 1.0; 0.5 |]
let of2 = [| 0; 10 |]
let dn6 = [| 0; 1; 3; 5; 7; 13 |]
let sc3 = [| 1.0; Float.ldexp 1.0 (-600); Float.ldexp 1.0 600 |]
let ofm = [| 0; 17 |]

let tgemm () =
  let ld = 8 and nb = 64 in
  let a = Array.make nb 0.0 and b = Array.make nb 0.0 and c = Array.make nb 0.0 in
  rseed 20261001;
  for ita = 1 to 3 do
    for itb = 1 to 3 do
      for im = 1 to 3 do
        for ix = 1 to 3 do
          for ik = 1 to 3 do
            for ia = 1 to 3 do
              for ib = 1 to 3 do
                for io = 1 to 2 do
                  let m = dm3.(im - 1) and n = dm3.(ix - 1) and k = dm3.(ik - 1) in
                  let o = of2.(io - 1) in
                  fillb a nb;
                  fillb b nb;
                  fillb c nb;
                  lbl 1 ita itb ((m * 100) + (n * 10) + k) ia ib io;
                  dgemm (tr ita) (tr itb) m n k al3.(ia - 1) a o ld b o ld
                    be3.(ib - 1) c o ld;
                  dumpb c nb
                done
              done
            done
          done
        done
      done
    done
  done

let ttri which =
  let ld = 8 and nb = 64 in
  let a = Array.make nb 0.0 and b = Array.make nb 0.0 in
  rseed (20261002 + which);
  for is = 1 to 2 do
    for iu = 1 to 2 do
      for it = 1 to 3 do
        for id = 1 to 2 do
          for im = 1 to 4 do
            for ix = 1 to 4 do
              for ia = 1 to 3 do
                for io = 1 to 2 do
                  let m = dm4.(im - 1) and n = dm4.(ix - 1) in
                  let o = of2.(io - 1) in
                  fillb a nb;
                  fillb b nb;
                  lbl (2 + which)
                    ((is * 1000) + (iu * 100) + (it * 10) + id)
                    ((m * 10) + n) ia io 0 0;
                  if which = 0 then
                    dtrsm (sd_ is) (up iu) (tr it) (dg id) m n al3.(ia - 1) a o
                      ld b o ld
                  else
                    dtrmm (sd_ is) (up iu) (tr it) (dg id) m n al3.(ia - 1) a o
                      ld b o ld;
                  dumpb b nb
                done
              done
            done
          done
        done
      done
    done
  done

let tsyrk () =
  let ld = 8 and nb = 64 in
  let a = Array.make nb 0.0 and c = Array.make nb 0.0 in
  rseed 20261004;
  for iu = 1 to 2 do
    for it = 1 to 3 do
      for ix = 1 to 4 do
        for ik = 1 to 4 do
          for ia = 1 to 3 do
            for ib = 1 to 3 do
              for io = 1 to 2 do
                let n = dm4.(ix - 1) and k = dm4.(ik - 1) in
                let o = of2.(io - 1) in
                fillb a nb;
                fillb c nb;
                lbl 4 ((iu * 10) + it) ((n * 10) + k) ia ib io 0;
                dsyrk (up iu) (tr it) n k al3.(ia - 1) a o ld be3.(ib - 1) c o
                  ld;
                dumpb c nb
              done
            done
          done
        done
      done
    done
  done

let tgemv () =
  let ld = 8 and nb = 64 in
  let a = Array.make nb 0.0 and x = Array.make nb 0.0 and y = Array.make nb 0.0 in
  let ic = [| 1; -2 |] in
  rseed 20261005;
  for it = 1 to 3 do
    for im = 1 to 4 do
      for ix = 1 to 4 do
        for jx = 1 to 2 do
          for jy = 1 to 2 do
            for ia = 1 to 3 do
              for ib = 1 to 3 do
                for io = 1 to 2 do
                  let m = dm4.(im - 1) and n = dm4.(ix - 1) in
                  let o = of2.(io - 1) in
                  fillb a nb;
                  fillb x nb;
                  fillb y nb;
                  lbl 5 it ((m * 10) + n) ((jx * 10) + jy) ia ib io;
                  dgemv (tr it) m n al3.(ia - 1) a o ld x o ic.(jx - 1)
                    be3.(ib - 1) y o ic.(jy - 1);
                  dumpb y nb
                done
              done
            done
          done
        done
      done
    done
  done

let tger () =
  let ld = 8 and nb = 64 in
  let a = Array.make nb 0.0 and x = Array.make nb 0.0 and y = Array.make nb 0.0 in
  let ic = [| 1; -2 |] in
  rseed 20261006;
  for im = 1 to 4 do
    for ix = 1 to 4 do
      for jx = 1 to 2 do
        for jy = 1 to 2 do
          for ia = 1 to 3 do
            for io = 1 to 2 do
              let m = dm4.(im - 1) and n = dm4.(ix - 1) in
              let o = of2.(io - 1) in
              fillb a nb;
              fillb x nb;
              fillb y nb;
              lbl 6 ((m * 10) + n) ((jx * 10) + jy) ia io 0 0;
              dger m n al3.(ia - 1) x o ic.(jx - 1) y o ic.(jy - 1) a o ld;
              dumpb a nb
            done
          done
        done
      done
    done
  done

let tvec () =
  let nb = 64 in
  let x = Array.make nb 0.0 and y = Array.make nb 0.0 in
  let ic = [| 1; 2; -1; -2 |] in
  rseed 20261007;
  for ix = 1 to 4 do
    for jx = 1 to 4 do
      for jy = 1 to 4 do
        for io = 1 to 2 do
          let n = dm4.(ix - 1) in
          let o = of2.(io - 1) in
          fillb x nb;
          fillb y nb;
          lbl 7 n ((jx * 10) + jy) io 0 0 0;
          dswap n x o ic.(jx - 1) y o ic.(jy - 1);
          dumpb x nb;
          dumpb y nb
        done
      done
    done
  done;
  for ix = 1 to 4 do
    for jx = 1 to 4 do
      for ia = 1 to 3 do
        for io = 1 to 2 do
          let n = dm4.(ix - 1) in
          let o = of2.(io - 1) in
          fillb x nb;
          lbl 8 n jx ia io 0 0;
          dscal n al3.(ia - 1) x o ic.(jx - 1);
          dumpb x nb
        done
      done
    done
  done;
  for ix = 1 to 4 do
    for jx = 1 to 4 do
      for io = 1 to 2 do
        let n = dm4.(ix - 1) in
        let o = of2.(io - 1) in
        fillb x nb;
        lbl 9 n jx io 0 0 0;
        pin (idamax n x o ic.(jx - 1) + 1)
      done
    done
  done;
  for ix = 1 to 6 do
    for jx = 1 to 4 do
      for jy = 1 to 4 do
        for io = 1 to 2 do
          let n = dn6.(ix - 1) in
          let o = of2.(io - 1) in
          fillb x nb;
          fillb y nb;
          lbl 50 n ((jx * 10) + jy) io 0 0 0;
          pr (ddot n x o ic.(jx - 1) y o ic.(jy - 1))
        done
      done
    done
  done;
  for ix = 1 to 6 do
    for jx = 1 to 4 do
      for io = 1 to 2 do
        let n = dn6.(ix - 1) in
        let o = of2.(io - 1) in
        fillb x nb;
        lbl 51 n jx io 0 0 0;
        pr (dasum n x o ic.(jx - 1))
      done
    done
  done;
  for ix = 1 to 6 do
    for jx = 1 to 4 do
      for is = 1 to 3 do
        for io = 1 to 2 do
          let n = dn6.(ix - 1) in
          let o = of2.(io - 1) in
          fillb x nb;
          let s = sc3.(is - 1) in
          let i = ref 1 in
          while !i < nb do
            x.(!i) <- x.(!i) *. s;
            i := !i + 2
          done;
          lbl 52 n ((jx * 10) + is) io 0 0 0;
          pr (dnrm2 n x o ic.(jx - 1))
        done
      done
    done
  done

let tlu which =
  let ld = 16 and nb = 256 in
  let a = Array.make nb 0.0 in
  let ipiv = Array.make ld 0 in
  let ms = [| 1; 1; 3; 2; 3; 5; 5; 3; 8; 13; 7 |] in
  let ns = [| 1; 3; 1; 3; 2; 5; 3; 5; 8; 7; 13 |] in
  rseed (20261010 + which);
  for ic = 1 to 11 do
    for io = 1 to 2 do
      for iz = 0 to 1 do
        let m = ms.(ic - 1) and n = ns.(ic - 1) in
        let o = ofm.(io - 1) in
        fillb a nb;
        if iz = 1 then begin
          let j = (n + 1) / 2 in
          for i = 1 to m do
            a.(o + i - 1 + ((j - 1) * ld)) <- 0.0
          done
        end;
        for i = 0 to ld - 1 do
          ipiv.(i) <- -8
        done;
        lbl (10 + which) m n io iz 0 0;
        let info =
          if which = 0 then dgetf2 m n a o ld ipiv 0
          else if which = 1 then dgetrf2 m n a o ld ipiv 0
          else dgetrf m n a o ld ipiv 0
        in
        pin info;
        for i = 0 to ld - 1 do
          pin (ipiv.(i) + 1)
        done;
        dumpb a nb
      done
    done
  done

let tchol which =
  let ld = 16 and nb = 256 in
  let a = Array.make nb 0.0 in
  let ns = [| 1; 2; 3; 5; 8; 13 |] in
  rseed (20261020 + which);
  for ic = 1 to 6 do
    for iu = 1 to 2 do
      for io = 1 to 2 do
        for iz = 0 to 1 do
          let n = ns.(ic - 1) in
          let o = ofm.(io - 1) in
          fillspd a nb o ld n;
          if iz = 1 then a.(o + (n - 1) + ((n - 1) * ld)) <- -1.0;
          lbl (20 + which) n iu io iz 0 0;
          let info =
            if which = 0 then dpotf2 (up iu) n a o ld
            else if which = 1 then dpotrf2 (up iu) n a o ld
            else dpotrf (up iu) n a o ld
          in
          pin info;
          dumpb a nb
        done
      done
    done
  done

let tsolve () =
  let ld = 16 and nb = 256 in
  let a = Array.make nb 0.0 and b = Array.make nb 0.0 in
  let ipiv = Array.make ld 0 in
  let ns = [| 1; 2; 3; 5; 8 |] in
  let rs = [| 1; 2; 3 |] in
  rseed 20261030;
  for ic = 1 to 5 do
    for ir = 1 to 3 do
      for io = 1 to 2 do
        let n = ns.(ic - 1) and nrhs = rs.(ir - 1) in
        let o = ofm.(io - 1) in
        fillb a nb;
        fillb b nb;
        for i = 0 to ld - 1 do
          ipiv.(i) <- -8
        done;
        lbl 30 n nrhs io 0 0 0;
        let info = dgesv n nrhs a o ld ipiv 0 b o ld in
        pin info;
        for i = 0 to ld - 1 do
          pin (ipiv.(i) + 1)
        done;
        dumpb a nb;
        dumpb b nb
      done
    done
  done;
  for ic = 1 to 5 do
    for ir = 1 to 3 do
      for it = 1 to 3 do
        for io = 1 to 2 do
          let n = ns.(ic - 1) and nrhs = rs.(ir - 1) in
          let o = ofm.(io - 1) in
          fillb a nb;
          fillb b nb;
          for i = 0 to ld - 1 do
            ipiv.(i) <- -8
          done;
          let fi = dgetrf n n a o ld ipiv 0 in
          lbl 31 n nrhs it io fi 0;
          let t = if it = 1 then No_trans else if it = 2 then Trans else Conj_trans in
          let info = dgetrs t n nrhs a o ld ipiv 0 b o ld in
          pin info;
          dumpb b nb
        done
      done
    done
  done;
  for ic = 1 to 5 do
    for ir = 1 to 3 do
      for iu = 1 to 2 do
        for io = 1 to 2 do
          let n = ns.(ic - 1) and nrhs = rs.(ir - 1) in
          let o = ofm.(io - 1) in
          fillspd a nb o ld n;
          fillb b nb;
          lbl 32 n nrhs iu io 0 0;
          let info = dposv (up iu) n nrhs a o ld b o ld in
          pin info;
          dumpb a nb;
          dumpb b nb
        done
      done
    done
  done;
  for ic = 1 to 5 do
    for ir = 1 to 3 do
      for iu = 1 to 2 do
        for io = 1 to 2 do
          let n = ns.(ic - 1) and nrhs = rs.(ir - 1) in
          let o = ofm.(io - 1) in
          fillspd a nb o ld n;
          fillb b nb;
          let fi = dpotrf (up iu) n a o ld in
          lbl 33 n nrhs iu io fi 0;
          let info = dpotrs (up iu) n nrhs a o ld b o ld in
          pin info;
          dumpb b nb
        done
      done
    done
  done

let tbig () =
  let ld = 104 and nb = 10816 in
  let a = Array.make nb 0.0 and b = Array.make nb 0.0 in
  let ipiv = Array.make ld 0 in
  let ns = [| 63; 64; 65; 70; 100 |] in
  rseed 20261040;
  for ic = 1 to 5 do
    let n = ns.(ic - 1) in
    fillb a nb;
    for i = 0 to ld - 1 do
      ipiv.(i) <- -8
    done;
    lbl 40 n 0 0 0 0 0;
    let info = dgetrf n n a 0 ld ipiv 0 in
    pin info;
    for i = 0 to ld - 1 do
      pin (ipiv.(i) + 1)
    done;
    dumpb a nb
  done;
  for ic = 1 to 5 do
    for iu = 1 to 2 do
      let n = ns.(ic - 1) in
      fillspd a nb 0 ld n;
      lbl 41 n iu 0 0 0 0;
      let info = dpotrf (up iu) n a 0 ld in
      pin info;
      dumpb a nb
    done
  done;
  for ic = 1 to 5 do
    let n = ns.(ic - 1) in
    fillb a nb;
    fillb b nb;
    for i = 0 to ld - 1 do
      ipiv.(i) <- -8
    done;
    lbl 42 n 0 0 0 0 0;
    let info = dgesv n 3 a 0 ld ipiv 0 b 0 ld in
    pin info;
    dumpb b nb
  done;
  for ic = 1 to 5 do
    for iu = 1 to 2 do
      let n = ns.(ic - 1) in
      fillspd a nb 0 ld n;
      fillb b nb;
      lbl 43 n iu 0 0 0 0;
      let info = dposv (up iu) n 3 a 0 ld b 0 ld in
      pin info;
      dumpb b nb
    done
  done

let () =
  tgemm ();
  ttri 0;
  ttri 1;
  tsyrk ();
  tgemv ();
  tger ();
  tvec ();
  tlu 0;
  tlu 1;
  tlu 2;
  tchol 0;
  tchol 1;
  tchol 2;
  tsolve ();
  tbig ();
  print_string (Buffer.contents buf)
