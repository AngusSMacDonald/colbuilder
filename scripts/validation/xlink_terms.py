"""Summarise the bonded terms that cross each crosslink bond in col_*.itp files.

Usage: python xlink_terms.py ITP_DIR [ITP_DIR2]

A crosslink bond is a [ bonds ] entry between two residues whose atom names
are not {C, N} (peptide bond). For each one, count the angles, proper
dihedrals, impropers and pairs that cross it, and record the (type, charge)
of each end. With two directories, a DIFF section lists the labels whose
distributions or end types/charges differ.
"""
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

SEC = re.compile(r"^\s*\[\s*(\w+)\s*\]")
PROPER = {1, 3, 9}
IMPROPER = {2, 4}
FIELDS = ("n", "ang", "prop", "imp", "pair", "ends")


def norm(a, b):
    return (a, b) if a < b else (b, a)


def parse_itp(path):
    atoms = {}  # nr -> (residue_key, resname, atomname, type, charge)
    bonds, pairs, angles, dihedrals = [], [], [], []
    section = None
    res_key, last_key = 0, None
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith("#"):
                continue
            m = SEC.match(line)
            if m:
                section = m.group(1)
                continue
            s = line.split(";", 1)[0].split()
            if not s:
                continue
            if section == "atoms":
                key = (s[2], s[3])
                if key != last_key:
                    res_key += 1
                    last_key = key
                atoms[int(s[0])] = (res_key, s[3], s[4], s[1], float(s[6]))
            elif section == "bonds":
                bonds.append((int(s[0]), int(s[1])))
            elif section == "pairs":
                pairs.append((int(s[0]), int(s[1])))
            elif section == "angles":
                angles.append((int(s[0]), int(s[1]), int(s[2])))
            elif section == "dihedrals":
                dihedrals.append((int(s[4]), (int(s[0]), int(s[1]), int(s[2]), int(s[3]))))
    return atoms, bonds, pairs, angles, dihedrals


def new_stat():
    return {"n": 0, "ang": Counter(), "prop": Counter(), "imp": Counter(),
            "pair": Counter(), "ends": (set(), set())}


def scan(path, total):
    """Add the crosslink stats of one itp into total; return its net charge."""
    atoms, bonds, pairs, angles, dihedrals = parse_itp(path)
    charge = sum(v[4] for v in atoms.values())

    adj = defaultdict(set)
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)

    xb = Counter()  # norm(a, b) -> number of [ bonds ] entries
    for a, b in bonds:
        ra, rb = atoms[a], atoms[b]
        if ra[0] == rb[0] or {ra[2], rb[2]} == {"C", "N"}:
            continue
        xb[norm(a, b)] += 1
    if not xb:
        return charge

    xb_keys = set(xb)
    atom_xb = defaultdict(set)  # atom nr -> crosslink keys containing it
    for a, b in xb:
        atom_xb[a].add((a, b))
        atom_xb[b].add((a, b))

    ang, prop, imp, prs = Counter(), Counter(), Counter(), Counter()

    # angle (i,j,k) crosses a bond only if its middle atom is a crosslink atom
    for i, j, k in angles:
        if j in atom_xb:
            for key in {norm(i, j), norm(j, k)} & xb_keys:
                ang[key] += 1

    for funct, quad in dihedrals:
        if not any(x in atom_xb for x in quad):
            continue
        i, j, k, l = quad
        if funct in PROPER:
            for key in {norm(i, j), norm(j, k), norm(k, l)} & xb_keys:
                prop[key] += 1
        elif funct in IMPROPER:
            if len({atoms[x][0] for x in quad}) < 2:
                continue
            hit = set()
            for x in quad:
                hit |= atom_xb.get(x, set())
            for key in hit:
                imp[key] += 1

    # pair (i,j) crosses a bond if a path i-x-y-j uses that bond as one edge;
    # such a path needs i or j within one bond of a crosslink atom
    near = set(atom_xb)
    for x in list(atom_xb):
        near.update(adj.get(x, ()))
    for i, j in pairs:
        if i not in near and j not in near:
            continue
        found = set()
        for x in adj.get(i, ()):
            if x == j:
                continue
            for y in adj.get(x, ()):
                if y in (i, j) or j not in adj.get(y, ()):
                    continue
                for e in (norm(i, x), norm(x, y), norm(y, j)):
                    if e in xb_keys:
                        found.add(e)
        for key in found:
            prs[key] += 1

    for key, n in xb.items():
        ia, ib = key
        ta = f"{atoms[ia][1]}:{atoms[ia][2]}"
        tb = f"{atoms[ib][1]}:{atoms[ib][2]}"
        label = "-".join(sorted([ta, tb]))
        st = total.setdefault(label, new_stat())
        ea = (atoms[ia][3], round(atoms[ia][4], 6))
        eb = (atoms[ib][3], round(atoms[ib][4], 6))
        for _ in range(n):
            st["n"] += 1
            st["ang"][ang[key]] += 1
            st["prop"][prop[key]] += 1
            st["imp"][imp[key]] += 1
            st["pair"][prs[key]] += 1
            if ta <= tb:
                st["ends"][0].add(ea)
                st["ends"][1].add(eb)
            else:
                st["ends"][0].add(eb)
                st["ends"][1].add(ea)
            if ta == tb:
                st["ends"][0].add(eb)
                st["ends"][1].add(ea)
    return charge


def fmt_counter(c):
    return "{" + ",".join(f"{k}:{v}" for k, v in sorted(c.items())) + "}"


def fmt_ends(s):
    return " ".join(f"{t}/{q:g}" for t, q in sorted(s))


def line(label, st):
    return (f"{label}  bonds={st['n']}  angles={fmt_counter(st['ang'])}"
            f" propers={fmt_counter(st['prop'])} impropers={fmt_counter(st['imp'])}"
            f" pairs={fmt_counter(st['pair'])}"
            f"  ends=[{fmt_ends(st['ends'][0])}] [{fmt_ends(st['ends'][1])}]")


def summarise(root):
    files = sorted(p for p in Path(root).glob("col_*.itp") if not p.name.startswith("posre"))
    total, charge = {}, 0.0
    for itp in files:
        charge += scan(itp, total)
    print(f"== {root}  ({len(files)} itp files)")
    print(f"net charge (sum of all itp charges): {charge:.4f}")
    print(f"crosslink bonds: {sum(st['n'] for st in total.values())}  labels: {len(total)}")
    for label in sorted(total):
        print(line(label, total[label]))
    print()
    return total


def diff(a, b, dir_a, dir_b):
    print("DIFF")
    print(f"A = {dir_a}")
    print(f"B = {dir_b}")
    count = 0
    for label in sorted(set(a) | set(b)):
        if label not in a or label not in b:
            count += 1
            print(f"{label}  only in {'A' if label in a else 'B'}")
            continue
        fields = [k for k in FIELDS if a[label][k] != b[label][k]]
        if fields:
            count += 1
            print(f"{label}  differs: {', '.join(fields)}")
            print("   A: " + line(label, a[label]))
            print("   B: " + line(label, b[label]))
    if not count:
        print("no differences")


def main(argv):
    if not 1 <= len(argv) <= 2:
        sys.exit(__doc__)
    results = [summarise(p) for p in argv]
    if len(results) == 2:
        diff(results[0], results[1], argv[0], argv[1])


if __name__ == "__main__":
    main(sys.argv[1:])
