# Summarises build/fort.txt against build/ocaml.txt: lines, how many were bit
# identical and, where they were not, the spread in ulps and the largest
# relative difference.  Naming a tag on the command line prints its first few
# differing lines.
#
# Near a zero the ulp count stops meaning much, because the result there is
# itself the difference of two nearly equal numbers, so the relative column is
# the one to read in that case.
import struct, sys, os, collections

here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(here, "build")


def load(p):
    with open(p) as f:
        return [(ln[:12].strip(), ln[12:20].strip(), ln[20:].rstrip("\n").strip())
                for ln in f]


def tofloat(h):
    return struct.unpack(">d", bytes.fromhex(h))[0]


def ordered(h):
    n = int(h, 16)
    return n if n < (1 << 63) else -(n - (1 << 63))


def pct(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))]


fa = load(os.path.join(out, "fort.txt"))
fb = load(os.path.join(out, "ocaml.txt"))
if len(fa) != len(fb):
    print("line counts differ: %d vs %d" % (len(fa), len(fb)))

tot = collections.Counter()
ulp = collections.defaultdict(list)
rel = collections.defaultdict(float)
examples = collections.defaultdict(list)
raised = collections.Counter()

for i, (ra, rb) in enumerate(zip(fa, fb)):
    tag = ra[0]
    tot[tag] += 1
    if ra == rb:
        continue
    if "RAISE" in ra[2] or "RAISE" in rb[2]:
        raised[tag] += 1
    if len(ra[2]) == 16 and len(rb[2]) == 16:
        x, y = tofloat(ra[2]), tofloat(rb[2])
        ulp[tag].append(abs(ordered(ra[2]) - ordered(rb[2])))
        d = max(abs(x), abs(y))
        if d not in (0.0, float("inf")) and d == d:
            rel[tag] = max(rel[tag], abs(x - y) / d)
    else:
        ulp[tag].append(-1)
    if len(examples[tag]) < 6:
        examples[tag].append((i + 1, ra[2], rb[2]))

print("%-12s %7s %7s %7s %7s %9s %11s" %
      ("tag", "lines", "differ", "med", "p90", "max ulp", "max rel"))
for tag in sorted(tot):
    u = ulp[tag]
    if not u:
        print("%-12s %7d %7d" % (tag, tot[tag], 0))
    else:
        print("%-12s %7d %7d %7d %7d %9d %11.3g" %
              (tag, tot[tag], len(u), pct(u, 0.5), pct(u, 0.9), max(u), rel[tag]))

print()
print("total lines %d, same %d, differing %d" %
      (len(fa), len(fa) - sum(len(v) for v in ulp.values()),
       sum(len(v) for v in ulp.values())))
if raised:
    print("RAISE on one side only:", dict(raised))

for tag in sys.argv[1:]:
    print()
    print("== %s ==" % tag)
    for ln, x, y in examples[tag]:
        print("line %d  fort %s  ocaml %s" % (ln, x, y))
