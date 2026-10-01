(lang dune 3.0)
(generate_opam_files true)

(name reeve)
(license Apache-2.0)
(authors "Zane Hambly")
(maintainers "Zane Hambly <zanehambly@gmail.com>")
(source (github Zaneham/reeve))

(package
 (name reeve)
 (synopsis "SLATEC and a few LAPACK routines in OCaml")
 (description
  "A port of the SLATEC Common Mathematical Library, with a few LAPACK routines, into OCaml. It is written in OCaml throughout, so the only thing it depends on is the compiler. It is checked by running it and the reference Fortran over the same inputs and comparing the output in hex.")
 (depends (ocaml (>= 4.14))))
