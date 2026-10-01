import re, sys, os

SRC = r"C:\dev\numerical\slatec-modern\src\modern\special_functions"
ORIG = r"C:\dev\numerical\slatec-modern\src\original\src"

def joincont(text):
    # join Fortran continuation lines (& at end)
    out = []
    buf = ""
    for line in text.splitlines():
        s = line.rstrip()
        if s.endswith("&"):
            buf += s[:-1]
        else:
            buf += s
            out.append(buf)
            buf = ""
    if buf:
        out.append(buf)
    return out

def find_array(fname, name):
    text = open(os.path.join(SRC, fname), encoding="utf-8", errors="replace").read()
    lines = joincont(text)
    pat = re.compile(r"\b" + re.escape(name) + r"\s*\(\s*[0-9,:]+\s*\)\s*=\s*(.*)$", re.I)
    for ln in lines:
        if ln.lstrip().startswith("!"):
            continue
        m = pat.search(ln)
        if m:
            rest = m.group(1)
            # take the bracketed list
            i = rest.find("[")
            if i < 0:
                continue
            depth = 0
            for j in range(i, len(rest)):
                if rest[j] == "[":
                    depth += 1
                elif rest[j] == "]":
                    depth -= 1
                    if depth == 0:
                        return rest[i+1:j]
            return rest[i+1:]
    return None

def to_ocaml(numstr):
    toks = [t.strip() for t in numstr.split(",")]
    toks = [t for t in toks if t]
    out = []
    for t in toks:
        t = re.sub(r"_DP\b", "", t, flags=re.I)
        t = t.strip()
        sign = ""
        if t.startswith("+"):
            t = t[1:]
        elif t.startswith("-"):
            sign = "-"
            t = t[1:]
        if t.startswith("."):
            t = "0" + t
        t = t.replace("E", "e").replace("D", "e")
        if "e" in t and "." not in t.split("e")[0]:
            a, b = t.split("e")
            t = a + ".0e" + b
        if "e" not in t and "." not in t:
            t = t + ".0"
        out.append(sign + t)
    return out

def emit(ocname, vals, perline=3):
    lines = []
    lines.append("let %s =" % ocname)
    lines.append("  [|")
    for i in range(0, len(vals), perline):
        chunk = vals[i:i+perline]
        lines.append("    " + "; ".join(chunk) + ";")
    # strip trailing ; on last
    lines[-1] = lines[-1].rstrip(";")
    lines.append("  |]")
    return "\n".join(lines) + "\n"

if __name__ == "__main__":
    spec = sys.argv[1:]
    for item in spec:
        fname, name, ocname, per = item.split(":")
        s = find_array(fname, name)
        if s is None:
            print("(* MISSING %s in %s *)" % (name, fname))
            continue
        vals = to_ocaml(s)
        print("(* %s %s n=%d *)" % (fname, name, len(vals)))
        print(emit(ocname, vals, int(per)))
