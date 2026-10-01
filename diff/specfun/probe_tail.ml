
let () =
  let p n v = Printf.printf "%-22s %.17g\n" n v in
  p "dai 0" (dai 0.0);
  p "dbi 0" (dbi 0.0);
  p "dai -1" (dai (-1.0));
  p "dai 1" (dai 1.0);
  p "dbi 1" (dbi 1.0);
  p "dai 5" (dai 5.0);
  p "dbi 5" (dbi 5.0);
  p "dai -10" (dai (-10.0));
  p "dbi -10" (dbi (-10.0));
  p "daie 10" (daie 10.0);
  p "dbie 10" (dbie 10.0);
  p "dai_xmax" dai_xmax;
  p "dbi_xmax" dbi_xmax;
  p "dai 100" (dai 100.0);
  p "dgamln 1.5" (dgamln 1.5);
  p "dgamln 2.5" (dgamln 2.5);
  p "dgamln 5" (dgamln 5.0);
  p "dgamln 0.5" (dgamln 0.5);
  p "dgamln 0.3" (dgamln 0.3);
  p "dgamln 4.5" (dgamln 4.5);
  p "dgamln 3.5" (dgamln 3.5);
  p "dgamln 101" (dgamln 101.0);
  p "dgamln 1000.5" (dgamln 1000.5);
  p "gamln_zmin" gamln_zmin;
  p "dpsi 1" (dpsi 1.0);
  p "dpsi 2" (dpsi 2.0);
  p "dpsi 2.5" (dpsi 2.5);
  p "dpsi 3.5" (dpsi 3.5);
  p "dpsi 0.25" (dpsi 0.25);
  p "dpsi 0.75" (dpsi 0.75);
  p "dpsi 20" (dpsi 20.0);
  p "dpsi -0.5" (dpsi (-0.5));
  p "dpsi -20.5" (dpsi (-20.5));
  p "drf 1 2 3" (drf 1.0 2.0 3.0);
  p "drf 0 1 1" (drf 0.0 1.0 1.0);
  p "drf 0 0.5 1" (drf 0.0 0.5 1.0);
  p "drf 0 1 2" (drf 0.0 1.0 2.0);
  p "drf 1 1 1" (drf 1.0 1.0 1.0);
  p "drc 1 1" (drc 1.0 1.0);
  p "drc 0 1" (drc 0.0 1.0);
  p "drc 1 2" (drc 1.0 2.0);
  p "drc 4 4" (drc 4.0 4.0);
  p "drd 1 1 1" (drd 1.0 1.0 1.0);
  p "drd 0 2 1" (drd 0.0 2.0 1.0);
  p "drd 1 2 3" (drd 1.0 2.0 3.0);
  p "drd 4 8 12" (drd 4.0 8.0 12.0);
  p "drf uplim corner" (drf 0.0 drf_lolim drf_uplim);
  p "drc uplim corner" (drc 0.0 drc_uplim);
  p "drd uplim corner" (drd 0.0 drd_lolim drd_uplim);
  p "drf errtol" drf_errtol;
  p "drd errtol" drd_errtol;
  p "drc errtol" drc_errtol;
  p "drd lolim" drd_lolim;
  p "drd uplim" drd_uplim;
  p "dcot pi/4" (dcot (pi /. 4.0));
  p "wronskian" (let h = 1e-6 in
     let aip = (dai (0.0 +. h) -. dai (0.0 -. h)) /. (2.0 *. h) in
     let bip = (dbi (0.0 +. h) -. dbi (0.0 -. h)) /. (2.0 *. h) in
     dai 0.0 *. bip -. aip *. dbi 0.0);
  p "1/pi" (1.0 /. pi)
