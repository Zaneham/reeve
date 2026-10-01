(* Reeve, Copyright 2026 Zane Hambly.
   BLAS level 2 and level 3 for dense matrices, ported from the reference
   Fortran. Matrices are column major with an explicit base offset and leading
   dimension, so the element (i, j) zero based lives at ao + i + j * lda. *)

type trans = No_trans | Trans | Conj_trans
type uplo = Upper | Lower
type side = Left | Right
type diag = Unit | Non_unit

let arg name pos = invalid_arg (name ^ ": argument " ^ string_of_int pos)
let is_notrans t = match t with No_trans -> true | Trans | Conj_trans -> false
let is_upper u = match u with Upper -> true | Lower -> false
let is_left s = match s with Left -> true | Right -> false
let is_nounit d = match d with Non_unit -> true | Unit -> false

let dgemv trans m n alpha a ao lda x xo incx beta y yo incy =
  if m < 0 then arg "dgemv" 2
  else if n < 0 then arg "dgemv" 3
  else if lda < max 1 m then arg "dgemv" 6
  else if incx = 0 then arg "dgemv" 8
  else if incy = 0 then arg "dgemv" 11
  else if m = 0 || n = 0 || (alpha = 0.0 && beta = 1.0) then ()
  else begin
    let lenx = if is_notrans trans then n else m in
    let leny = if is_notrans trans then m else n in
    let kx = if incx > 0 then 0 else (1 - lenx) * incx in
    let ky = if incy > 0 then 0 else (1 - leny) * incy in
    let aij i j = a.(ao + i - 1 + ((j - 1) * lda)) in
    let xi i = x.(xo + kx + ((i - 1) * incx)) in
    let yi i = y.(yo + ky + ((i - 1) * incy)) in
    let sety i v = y.(yo + ky + ((i - 1) * incy)) <- v in
    if beta <> 1.0 then
      if beta = 0.0 then
        for i = 1 to leny do
          sety i 0.0
        done
      else
        for i = 1 to leny do
          sety i (beta *. yi i)
        done;
    if alpha <> 0.0 then
      if is_notrans trans then
        for j = 1 to n do
          let t = alpha *. xi j in
          for i = 1 to m do
            sety i (yi i +. (t *. aij i j))
          done
        done
      else
        for j = 1 to n do
          let t = ref 0.0 in
          for i = 1 to m do
            t := !t +. (aij i j *. xi i)
          done;
          sety j (yi j +. (alpha *. !t))
        done
  end

let dger m n alpha x xo incx y yo incy a ao lda =
  if m < 0 then arg "dger" 1
  else if n < 0 then arg "dger" 2
  else if incx = 0 then arg "dger" 5
  else if incy = 0 then arg "dger" 7
  else if lda < max 1 m then arg "dger" 9
  else if m = 0 || n = 0 || alpha = 0.0 then ()
  else begin
    let jy = if incy > 0 then 0 else (1 - n) * incy in
    let kx = if incx > 0 then 0 else (1 - m) * incx in
    for j = 1 to n do
      let yj = y.(yo + jy + ((j - 1) * incy)) in
      if yj <> 0.0 then begin
        let t = alpha *. yj in
        for i = 1 to m do
          let p = ao + i - 1 + ((j - 1) * lda) in
          a.(p) <- a.(p) +. (x.(xo + kx + ((i - 1) * incx)) *. t)
        done
      end
    done
  end

let dgemm transa transb m n k alpha a ao lda b bo ldb beta c co ldc =
  let nota = is_notrans transa in
  let notb = is_notrans transb in
  let nrowa = if nota then m else k in
  let nrowb = if notb then k else n in
  if m < 0 then arg "dgemm" 3
  else if n < 0 then arg "dgemm" 4
  else if k < 0 then arg "dgemm" 5
  else if lda < max 1 nrowa then arg "dgemm" 8
  else if ldb < max 1 nrowb then arg "dgemm" 10
  else if ldc < max 1 m then arg "dgemm" 13
  else if m = 0 || n = 0 || ((alpha = 0.0 || k = 0) && beta = 1.0) then ()
  else begin
    let aij i j = a.(ao + i - 1 + ((j - 1) * lda)) in
    let bij i j = b.(bo + i - 1 + ((j - 1) * ldb)) in
    let cij i j = c.(co + i - 1 + ((j - 1) * ldc)) in
    let setc i j v = c.(co + i - 1 + ((j - 1) * ldc)) <- v in
    if alpha = 0.0 then
      if beta = 0.0 then
        for j = 1 to n do
          for i = 1 to m do
            setc i j 0.0
          done
        done
      else
        for j = 1 to n do
          for i = 1 to m do
            setc i j (beta *. cij i j)
          done
        done
    else if notb then
      if nota then
        for j = 1 to n do
          if beta = 0.0 then
            for i = 1 to m do
              setc i j 0.0
            done
          else if beta <> 1.0 then
            for i = 1 to m do
              setc i j (beta *. cij i j)
            done;
          for l = 1 to k do
            let t = alpha *. bij l j in
            for i = 1 to m do
              setc i j (cij i j +. (t *. aij i l))
            done
          done
        done
      else
        for j = 1 to n do
          for i = 1 to m do
            let t = ref 0.0 in
            for l = 1 to k do
              t := !t +. (aij l i *. bij l j)
            done;
            if beta = 0.0 then setc i j (alpha *. !t)
            else setc i j ((alpha *. !t) +. (beta *. cij i j))
          done
        done
    else if nota then
      for j = 1 to n do
        if beta = 0.0 then
          for i = 1 to m do
            setc i j 0.0
          done
        else if beta <> 1.0 then
          for i = 1 to m do
            setc i j (beta *. cij i j)
          done;
        for l = 1 to k do
          let t = alpha *. bij j l in
          for i = 1 to m do
            setc i j (cij i j +. (t *. aij i l))
          done
        done
      done
    else
      for j = 1 to n do
        for i = 1 to m do
          let t = ref 0.0 in
          for l = 1 to k do
            t := !t +. (aij l i *. bij j l)
          done;
          if beta = 0.0 then setc i j (alpha *. !t)
          else setc i j ((alpha *. !t) +. (beta *. cij i j))
        done
      done
  end

let dsyrk uplo trans n k alpha a ao lda beta c co ldc =
  let notr = is_notrans trans in
  let upper = is_upper uplo in
  let nrowa = if notr then n else k in
  if n < 0 then arg "dsyrk" 3
  else if k < 0 then arg "dsyrk" 4
  else if lda < max 1 nrowa then arg "dsyrk" 7
  else if ldc < max 1 n then arg "dsyrk" 10
  else if n = 0 || ((alpha = 0.0 || k = 0) && beta = 1.0) then ()
  else begin
    let aij i j = a.(ao + i - 1 + ((j - 1) * lda)) in
    let cij i j = c.(co + i - 1 + ((j - 1) * ldc)) in
    let setc i j v = c.(co + i - 1 + ((j - 1) * ldc)) <- v in
    let lo j = if upper then 1 else j in
    let hi j = if upper then j else n in
    if alpha = 0.0 then
      if beta = 0.0 then
        for j = 1 to n do
          for i = lo j to hi j do
            setc i j 0.0
          done
        done
      else
        for j = 1 to n do
          for i = lo j to hi j do
            setc i j (beta *. cij i j)
          done
        done
    else if notr then
      for j = 1 to n do
        if beta = 0.0 then
          for i = lo j to hi j do
            setc i j 0.0
          done
        else if beta <> 1.0 then
          for i = lo j to hi j do
            setc i j (beta *. cij i j)
          done;
        for l = 1 to k do
          if aij j l <> 0.0 then begin
            let t = alpha *. aij j l in
            for i = lo j to hi j do
              setc i j (cij i j +. (t *. aij i l))
            done
          end
        done
      done
    else
      for j = 1 to n do
        for i = lo j to hi j do
          let t = ref 0.0 in
          for l = 1 to k do
            t := !t +. (aij l i *. aij l j)
          done;
          if beta = 0.0 then setc i j (alpha *. !t)
          else setc i j ((alpha *. !t) +. (beta *. cij i j))
        done
      done
  end

let dtrsm side uplo transa diag m n alpha a ao lda b bo ldb =
  let left = is_left side in
  let upper = is_upper uplo in
  let notr = is_notrans transa in
  let nounit = is_nounit diag in
  let nrowa = if left then m else n in
  if m < 0 then arg "dtrsm" 5
  else if n < 0 then arg "dtrsm" 6
  else if lda < max 1 nrowa then arg "dtrsm" 9
  else if ldb < max 1 m then arg "dtrsm" 11
  else if m = 0 || n = 0 then ()
  else begin
    let aij i j = a.(ao + i - 1 + ((j - 1) * lda)) in
    let bij i j = b.(bo + i - 1 + ((j - 1) * ldb)) in
    let setb i j v = b.(bo + i - 1 + ((j - 1) * ldb)) <- v in
    if alpha = 0.0 then
      for j = 1 to n do
        for i = 1 to m do
          setb i j 0.0
        done
      done
    else if left then
      if notr then
        if upper then
          for j = 1 to n do
            if alpha <> 1.0 then
              for i = 1 to m do
                setb i j (alpha *. bij i j)
              done;
            for k = m downto 1 do
              if bij k j <> 0.0 then begin
                if nounit then setb k j (bij k j /. aij k k);
                for i = 1 to k - 1 do
                  setb i j (bij i j -. (bij k j *. aij i k))
                done
              end
            done
          done
        else
          for j = 1 to n do
            if alpha <> 1.0 then
              for i = 1 to m do
                setb i j (alpha *. bij i j)
              done;
            for k = 1 to m do
              if bij k j <> 0.0 then begin
                if nounit then setb k j (bij k j /. aij k k);
                for i = k + 1 to m do
                  setb i j (bij i j -. (bij k j *. aij i k))
                done
              end
            done
          done
      else if upper then
        for j = 1 to n do
          for i = 1 to m do
            let t = ref (alpha *. bij i j) in
            for k = 1 to i - 1 do
              t := !t -. (aij k i *. bij k j)
            done;
            if nounit then t := !t /. aij i i;
            setb i j !t
          done
        done
      else
        for j = 1 to n do
          for i = m downto 1 do
            let t = ref (alpha *. bij i j) in
            for k = i + 1 to m do
              t := !t -. (aij k i *. bij k j)
            done;
            if nounit then t := !t /. aij i i;
            setb i j !t
          done
        done
    else if notr then
      if upper then
        for j = 1 to n do
          if alpha <> 1.0 then
            for i = 1 to m do
              setb i j (alpha *. bij i j)
            done;
          for k = 1 to j - 1 do
            if aij k j <> 0.0 then
              for i = 1 to m do
                setb i j (bij i j -. (aij k j *. bij i k))
              done
          done;
          if nounit then begin
            let t = 1.0 /. aij j j in
            for i = 1 to m do
              setb i j (t *. bij i j)
            done
          end
        done
      else
        for j = n downto 1 do
          if alpha <> 1.0 then
            for i = 1 to m do
              setb i j (alpha *. bij i j)
            done;
          for k = j + 1 to n do
            if aij k j <> 0.0 then
              for i = 1 to m do
                setb i j (bij i j -. (aij k j *. bij i k))
              done
          done;
          if nounit then begin
            let t = 1.0 /. aij j j in
            for i = 1 to m do
              setb i j (t *. bij i j)
            done
          end
        done
    else if upper then
      for k = n downto 1 do
        if nounit then begin
          let t = 1.0 /. aij k k in
          for i = 1 to m do
            setb i k (t *. bij i k)
          done
        end;
        for j = 1 to k - 1 do
          if aij j k <> 0.0 then begin
            let t = aij j k in
            for i = 1 to m do
              setb i j (bij i j -. (t *. bij i k))
            done
          end
        done;
        if alpha <> 1.0 then
          for i = 1 to m do
            setb i k (alpha *. bij i k)
          done
      done
    else
      for k = 1 to n do
        if nounit then begin
          let t = 1.0 /. aij k k in
          for i = 1 to m do
            setb i k (t *. bij i k)
          done
        end;
        for j = k + 1 to n do
          if aij j k <> 0.0 then begin
            let t = aij j k in
            for i = 1 to m do
              setb i j (bij i j -. (t *. bij i k))
            done
          end
        done;
        if alpha <> 1.0 then
          for i = 1 to m do
            setb i k (alpha *. bij i k)
          done
      done
  end

let dtrmm side uplo transa diag m n alpha a ao lda b bo ldb =
  let left = is_left side in
  let upper = is_upper uplo in
  let notr = is_notrans transa in
  let nounit = is_nounit diag in
  let nrowa = if left then m else n in
  if m < 0 then arg "dtrmm" 5
  else if n < 0 then arg "dtrmm" 6
  else if lda < max 1 nrowa then arg "dtrmm" 9
  else if ldb < max 1 m then arg "dtrmm" 11
  else if m = 0 || n = 0 then ()
  else begin
    let aij i j = a.(ao + i - 1 + ((j - 1) * lda)) in
    let bij i j = b.(bo + i - 1 + ((j - 1) * ldb)) in
    let setb i j v = b.(bo + i - 1 + ((j - 1) * ldb)) <- v in
    if alpha = 0.0 then
      for j = 1 to n do
        for i = 1 to m do
          setb i j 0.0
        done
      done
    else if left then
      if notr then
        if upper then
          for j = 1 to n do
            for k = 1 to m do
              if bij k j <> 0.0 then begin
                let t = ref (alpha *. bij k j) in
                for i = 1 to k - 1 do
                  setb i j (bij i j +. (!t *. aij i k))
                done;
                if nounit then t := !t *. aij k k;
                setb k j !t
              end
            done
          done
        else
          for j = 1 to n do
            for k = m downto 1 do
              if bij k j <> 0.0 then begin
                let t = alpha *. bij k j in
                setb k j t;
                if nounit then setb k j (bij k j *. aij k k);
                for i = k + 1 to m do
                  setb i j (bij i j +. (t *. aij i k))
                done
              end
            done
          done
      else if upper then
        for j = 1 to n do
          for i = m downto 1 do
            let t = ref (bij i j) in
            if nounit then t := !t *. aij i i;
            for k = 1 to i - 1 do
              t := !t +. (aij k i *. bij k j)
            done;
            setb i j (alpha *. !t)
          done
        done
      else
        for j = 1 to n do
          for i = 1 to m do
            let t = ref (bij i j) in
            if nounit then t := !t *. aij i i;
            for k = i + 1 to m do
              t := !t +. (aij k i *. bij k j)
            done;
            setb i j (alpha *. !t)
          done
        done
    else if notr then
      if upper then
        for j = n downto 1 do
          let t = if nounit then alpha *. aij j j else alpha in
          for i = 1 to m do
            setb i j (t *. bij i j)
          done;
          for k = 1 to j - 1 do
            if aij k j <> 0.0 then begin
              let u = alpha *. aij k j in
              for i = 1 to m do
                setb i j (bij i j +. (u *. bij i k))
              done
            end
          done
        done
      else
        for j = 1 to n do
          let t = if nounit then alpha *. aij j j else alpha in
          for i = 1 to m do
            setb i j (t *. bij i j)
          done;
          for k = j + 1 to n do
            if aij k j <> 0.0 then begin
              let u = alpha *. aij k j in
              for i = 1 to m do
                setb i j (bij i j +. (u *. bij i k))
              done
            end
          done
        done
    else if upper then
      for k = 1 to n do
        for j = 1 to k - 1 do
          if aij j k <> 0.0 then begin
            let u = alpha *. aij j k in
            for i = 1 to m do
              setb i j (bij i j +. (u *. bij i k))
            done
          end
        done;
        let t = if nounit then alpha *. aij k k else alpha in
        if t <> 1.0 then
          for i = 1 to m do
            setb i k (t *. bij i k)
          done
      done
    else
      for k = n downto 1 do
        for j = k + 1 to n do
          if aij j k <> 0.0 then begin
            let u = alpha *. aij j k in
            for i = 1 to m do
              setb i j (bij i j +. (u *. bij i k))
            done
          end
        done;
        let t = if nounit then alpha *. aij k k else alpha in
        if t <> 1.0 then
          for i = 1 to m do
            setb i k (t *. bij i k)
          done
      done
  end
