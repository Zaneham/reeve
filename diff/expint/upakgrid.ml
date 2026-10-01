let hex v = Printf.sprintf "%016LX" (Int64.bits_of_float v)
let inp = open_out "upak.in"
let out = open_out "upak.ml"
let emit x =
  Printf.fprintf inp "d9upak %s 0\n" (hex x);
  let y, n = Float.frexp x in
  Printf.fprintf out "d9upak_n %s %d\nd9upak %s\n" (hex y) n (hex y)
let () =
  let st = Random.State.make [| 20261001 |] in
  emit 0.0; emit (-0.0);
  for e = -1074 to 1023 do
    emit (Float.ldexp 1.0 e);
    emit (-.Float.ldexp 1.0 e);
    emit (Float.ldexp 1.9999999999999998 e)
  done;
  let i = ref 0 in
  while !i < 200000 do
    let b = Random.State.int64 st Int64.max_int in
    let v = Int64.float_of_bits b in
    if Float.is_finite v then begin emit v; emit (-.v); incr i end
  done;
  close_out inp; close_out out
