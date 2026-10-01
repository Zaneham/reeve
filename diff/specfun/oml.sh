#!/bin/sh
# oml.sh <specfun.ml> <name>  -> one value per line
awk -v n="$2" '
  $0 == "let " n " =" { inside=1; next }
  inside && /^  \|\]/ { exit }
  inside && /;/ { gsub(/[ ;]/,""); if ($0 != "[|" && $0 != "") print }
' "$1"
