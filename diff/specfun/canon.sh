#!/bin/sh
f="$1"; name="$2"
tr -d '\r' < "$f" | tr -d '\n' \
  | sed "s/.*:: *${name}([0-9]*) *= *\[//" \
  | sed 's/\].*//' \
  | sed 's/&//g; s/_DP//g; s/_SP//g; s/ //g' \
  | tr ',' '\n' | grep -v '^$' \
  | awk '{ s=$0; sub(/^\+/,"",s); gsub(/E/,"e",s);
           if (substr(s,1,1)==".") s="0" s;
           else if (substr(s,1,2)=="-.") s="-0" substr(s,2);
           print s }'
