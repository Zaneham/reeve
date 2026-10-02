(* Reeve, Copyright 2026 Zane Hambly.
   The dense matrix the library hands around. *)

type t = { data : float array; rows : int; cols : int }

let create rows cols =
  if rows < 0 || cols < 0 then invalid_arg "Mat.create: negative dimension";
  { data = Array.make (max 1 (rows * cols)) 0.0; rows; cols }

let of_array rows cols data =
  if rows < 0 || cols < 0 then invalid_arg "Mat.of_array: negative dimension";
  if Array.length data < rows * cols then
    invalid_arg "Mat.of_array: array shorter than rows * cols";
  { data; rows; cols }

let of_vec v = { data = v; rows = Array.length v; cols = 1 }
let rows m = m.rows
let cols m = m.cols
let data m = m.data
let get m i j = m.data.(i + (j * m.rows))
let set m i j v = m.data.(i + (j * m.rows)) <- v
let copy m = { m with data = Array.copy m.data }

let init rows cols f =
  let m = create rows cols in
  for j = 0 to cols - 1 do
    for i = 0 to rows - 1 do
      set m i j (f i j)
    done
  done;
  m

let square name m =
  if m.rows <> m.cols then invalid_arg (name ^ ": matrix is not square")

let same_rows name a b =
  if a.rows <> b.rows then
    invalid_arg (name ^ ": operands disagree on the number of rows")
