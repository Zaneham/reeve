# Pairs the two sweeps case by case and says what every difference is. It
# reads work/in.*, work/f.out and work/ml.out from the directory it is run in.

import collections
import glob
import sys


def take(lines, i):
    head = lines[i]
    n = 0
    if head.startswith(("legf ", "nrmp ")) and not head.endswith("err"):
        n = int(head.split()[1])
    return lines[i : i + 1 + n], i + 1 + n


def negzero_only(bf, bm):
    """True when the two blocks differ only in the sign of a zero."""
    if len(bf) != len(bm):
        return False
    hit = False
    for a, b in zip(bf, bm):
        if a == b:
            continue
        pa, pb = a.split(), b.split()
        if len(pa) != len(pb) or pa[0] != pb[0]:
            return False
        for u, v in zip(pa, pb):
            if u == v:
                continue
            if {u, v} == {"8000000000000000", "0000000000000000"}:
                hit = True
            else:
                return False
    return hit


def main():
    ml = open("ml.out").read().split("\n")
    f = open("f.out").read().split("\n")
    im = iff = 0
    total = reachable = 0
    cause = collections.Counter()
    first = {}
    desync = None
    for name in sorted(glob.glob("in.*")):
        cmds = [c for c in open(name).read().split("\n")[1:] if c.strip()]
        total += len(cmds)
        if ml[im] != f[iff]:
            print(f"{name}: dxset disagrees, {f[iff]!r} against {ml[im]!r}")
            cause["dxset"] += 1
        ok = ml[im] == "set ok"
        im += 1
        iff += 1
        if not ok:
            continue
        reachable += len(cmds)
        for c in cmds:
            bm, im = take(ml, im)
            bf, iff = take(f, iff)
            if bm == bf:
                continue
            if len(bm) != len(bf):
                k = "length mismatch, the streams are now out of step"
                if desync is None:
                    desync = (name, c, bf, bm)
            elif negzero_only(bf, bm):
                k = "the sign of a zero"
            else:
                k = "the value"
            cause[k] += 1
            first.setdefault(k, (name, c, bf, bm))
    print(f"cases written {total}, reachable {reachable}")
    if not cause:
        print("no differences")
        return 0
    for k, v in cause.most_common():
        print(f"{v:8d}  {k}")
    for k, (name, c, bf, bm) in first.items():
        print(f"\nfirst '{k}' difference, {name}: {c}")
        for a, b in zip(bf, bm):
            mark = "  " if a == b else "> "
            print(f"  {mark}F  {a}")
            print(f"  {mark}ML {b}")
    return 1


sys.exit(main())
