"""Check local geometry around collagen crosslinks between two .gro files of one system.

usage: python xlink_strain.py TOP_FILE ITP_DIR FF_DIR GRO_BEFORE GRO_AFTER

Check A: bond-angle strain of [ angles ] terms touching a crosslink bond, against
         ffbonded.itp [ angletypes ] theta0.
Check B: out-of-plane distance of the pyridinium ring atoms (LYX/LXX/LXY/LYY).
Check C: chirality flips (signed volume of centre with 4 bonded neighbours) between
         the two structures, plus a CA handedness count in GRO_AFTER. Centres with
         two or more hydrogen neighbours (CH2/CH3) are skipped, since their
         hydrogens are interchangeable and the sign is not defined.

Global atom indices are 0-based, in [ molecules ] order times count.
Coordinates are converted from nm to Angstrom.
"""
import os
import re
import sys
from collections import Counter, defaultdict
from itertools import islice

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unbonded_xlinks import parse_top_molecules  # noqa: E402

SECTION_RE = re.compile(r"^\s*\[\s*(\S+?)\s*\]")

# Crosslink atom-name pairs; the two atoms must also sit in different residues
XLINK_NAME_PAIRS = {
    frozenset(("CE", "NZ")),
    frozenset(("C12", "CB")),
    frozenset(("C12", "CG")),
    frozenset(("C13", "CG")),
    frozenset(("CE", "NZ1")),
    frozenset(("CE1",)),
}
RING_RES = {"LYX", "LXX", "LXY", "LYY"}
RING_ATOMS = ("NZ", "C10", "C11", "C12", "C13", "C14")
CHIRAL_RES = {
    "LYX", "LXX", "LXY", "LYY", "LY2", "LY3", "LX2", "LX3",
    "L2Y", "L3Y", "L2X", "L3X", "L4Y", "L5Y", "L4X", "L5X",
    "LY4", "LY5", "LX4", "LX5",
}
MAX_EXAMPLES = 10
TOP_WORST = 5


def parse_itp_full(itp_file):
    """Return (atoms, bonds, angles) for one ITP file.

    atoms: list of (type, resnr, resname, atomname), index = local nr - 1
    bonds, angles: lists of 1-based local atom numbers
    Sections may repeat; everything is accumulated.
    """
    atoms_by_nr = {}
    bonds, angles = [], []
    section = None
    with open(itp_file) as fh:
        for raw in fh:
            line = raw.split(";", 1)[0].strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("["):
                m = SECTION_RE.match(line)
                section = m.group(1) if m else None
                continue
            parts = line.split()
            if section == "atoms":
                # nr type resnr residue atom cgnr charge mass
                atoms_by_nr[int(parts[0])] = (parts[1], int(parts[2]), parts[3], parts[4])
            elif section == "bonds":
                bonds.append((int(parts[0]), int(parts[1])))
            elif section == "angles":
                angles.append((int(parts[0]), int(parts[1]), int(parts[2])))
    atoms = [atoms_by_nr[nr] for nr in sorted(atoms_by_nr)]
    return atoms, bonds, angles


def molecule_info(itp_file):
    """Parse one ITP and derive its crosslink bonds, crosslink angles and rings."""
    atoms, bonds, angles = parse_itp_full(itp_file)

    xl_bonds = set()
    for a, b in bonds:
        ta, tb = atoms[a - 1], atoms[b - 1]
        if ta[2] != tb[2] and frozenset((ta[3], tb[3])) in XLINK_NAME_PAIRS:
            xl_bonds.add(frozenset((a, b)))
    xl_angles = [
        (i, j, k) for i, j, k in angles
        if frozenset((i, j)) in xl_bonds or frozenset((j, k)) in xl_bonds
    ]

    ring_atoms = defaultdict(dict)  # resnr -> {atomname: 0-based local index}
    for idx, (_, resnr, resname, name) in enumerate(atoms):
        if resname in RING_RES:
            ring_atoms[resnr][name] = idx
    rings = []
    for names in ring_atoms.values():
        if all(name in names for name in RING_ATOMS):
            rings.append([names[name] for name in RING_ATOMS])
    incomplete = len(ring_atoms) - len(rings)

    return {
        "atoms": atoms,
        "bonds": bonds,
        "xl_angles": xl_angles,
        "n_xl_bonds": len(xl_bonds),
        "rings": rings,
        "incomplete": incomplete,
    }


def build_topology(top_file, itp_dir):
    """Return (cache, instances, names, resnames, types, mols) in global order.

    instances: list of (offset, molname), one per molecule copy.
    """
    cache = {}
    names, resnames, types, mols = [], [], [], []
    instances = []
    for molname, count in parse_top_molecules(top_file):
        itp_file = os.path.join(itp_dir, molname + ".itp")
        if not os.path.isfile(itp_file):
            sys.exit(f"ERROR: no ITP file for molecule {molname}: {itp_file}")
        if molname not in cache:
            cache[molname] = molecule_info(itp_file)
        info = cache[molname]
        for _ in range(count):
            instances.append((len(names), molname))
            for typ, _, resname, name in info["atoms"]:
                types.append(typ)
                resnames.append(resname)
                names.append(name)
                mols.append(molname)
    return cache, instances, names, resnames, types, mols


def load_angletypes(ff_dir):
    """Return {(type_i, type_j, type_k): th0} from ffbonded.itp [ angletypes ]."""
    table = {}
    section = None
    with open(os.path.join(ff_dir, "ffbonded.itp")) as fh:
        for raw in fh:
            line = raw.split(";", 1)[0].strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("["):
                m = SECTION_RE.match(line)
                section = m.group(1) if m else None
                continue
            if section == "angletypes":
                parts = line.split()
                # i j k func th0 cth
                table.setdefault(tuple(parts[0:3]), float(parts[4]))
    return table


def theta0_for(table, ti, tj, tk):
    th0 = table.get((ti, tj, tk))
    if th0 is None:
        th0 = table.get((tk, tj, ti))
    return th0


def load_gro(gro_file, resnames, names):
    """Return xyz (Angstrom) for the first len(names) atoms, after checking names match."""
    n = len(names)
    with open(gro_file) as fh:
        fh.readline()  # title line
        n_gro = int(fh.readline().strip())
        if n_gro < n:
            sys.exit(f"SANITY FAIL: {gro_file} has {n_gro} atoms, topology needs {n}")
        gro_res, gro_name, xs, ys, zs = [], [], [], [], []
        lines = list(islice(fh, n))
    if len(lines) != n:
        sys.exit(f"SANITY FAIL: {gro_file} ended after {len(lines)} atom lines")
    for line in lines:
        gro_res.append(line[5:10].strip())
        gro_name.append(line[10:15].strip())
        xs.append(float(line[20:28]))
        ys.append(float(line[28:36]))
        zs.append(float(line[36:44]))
    if gro_res != resnames or gro_name != names:
        for i in range(n):
            if gro_res[i] != resnames[i] or gro_name[i] != names[i]:
                sys.exit(
                    f"SANITY FAIL: {gro_file} atom #{i + 1} is {gro_res[i]}:{gro_name[i]}, "
                    f"topology says {resnames[i]}:{names[i]}"
                )
    return np.column_stack((xs, ys, zs)) * 10.0


def angles_deg(xyz, idx):
    a = xyz[idx[:, 0]] - xyz[idx[:, 1]]
    b = xyz[idx[:, 2]] - xyz[idx[:, 1]]
    cos = np.einsum("ij,ij->i", a, b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
    return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))


def crosslink_angles(cache, instances):
    """Return (m, 3) global 0-based angle triples that touch a crosslink bond."""
    rows = []
    for offset, molname in instances:
        for i, j, k in cache[molname]["xl_angles"]:
            rows.append((offset + i - 1, offset + j - 1, offset + k - 1))
    return np.array(rows, dtype=int).reshape(-1, 3)


def check_angles(label, gro_file, xyz, idx, theta0, mols, resnames, names, print_worst):
    ang = angles_deg(xyz, idx)
    param = ~np.isnan(theta0)
    dev_all = np.abs(ang - theta0)
    dev = dev_all[param]
    print(f"{label}: {gro_file}")
    print(f"  crosslink angles={len(idx)} parameterised={param.sum()} "
          f"unparameterised={(~param).sum()}")
    if dev.size:
        print(f"  mean |dev|={dev.mean():.2f} max |dev|={dev.max():.2f} "
              f"|dev|>15: {(dev > 15).sum()} |dev|>25: {(dev > 25).sum()}")
    if print_worst:
        order = np.argsort(-np.where(param, dev_all, -1.0))
        print(f"  {TOP_WORST} worst (parameterised):")
        for r in order[:TOP_WORST]:
            if not param[r]:
                break
            i, j, k = idx[r]
            label_t = (f"{resnames[i]}:{names[i]}-{resnames[j]}:{names[j]}-"
                       f"{resnames[k]}:{names[k]}")
            print(f"    {mols[i]} {label_t} {ang[r]:.2f} {theta0[r]:.2f}")


def ring_oop(xyz, rings):
    """Max absolute out-of-plane distance (Angstrom) per ring, from an SVD plane fit."""
    pts = xyz[rings]                                  # (R, 6, 3)
    centred = pts - pts.mean(axis=1, keepdims=True)
    _, _, vt = np.linalg.svd(centred)
    normal = vt[:, 2, :]                              # smallest singular vector
    dist = np.einsum("rij,rj->ri", centred, normal)
    return np.abs(dist).max(axis=1)


def signed_volumes(xyz, centre, nb):
    """v = dot(n1-c, cross(n2-c, n3-c)) for each row; nb is (m, 3)."""
    c = xyz[centre]
    a = xyz[nb[:, 0]] - c
    b = xyz[nb[:, 1]] - c
    d = xyz[nb[:, 2]] - c
    return np.einsum("ij,ij->i", a, np.cross(b, d))


def main():
    if len(sys.argv) != 6:
        print("usage: python xlink_strain.py TOP_FILE ITP_DIR FF_DIR GRO_BEFORE GRO_AFTER",
              file=sys.stderr)
        sys.exit(2)
    top_file, itp_dir, ff_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    gro_before, gro_after = sys.argv[4], sys.argv[5]

    cache, instances, names, resnames, types, mols = build_topology(top_file, itp_dir)
    n_top = len(names)
    n_xl = sum(cache[m]["n_xl_bonds"] for _, m in instances)
    print(f"topology: atoms={n_top} molecule_copies={len(instances)} crosslink_bonds={n_xl}")

    xyz_b = load_gro(gro_before, resnames, names)
    xyz_a = load_gro(gro_after, resnames, names)

    # ---- Check A: angle strain ---------------------------------------
    print("\n== Check A: crosslink angle strain")
    table = load_angletypes(ff_dir)
    idx = crosslink_angles(cache, instances)
    theta0 = np.array([
        np.nan if (t := theta0_for(table, types[i], types[j], types[k])) is None else t
        for i, j, k in idx
    ], dtype=float)
    unparam = Counter(
        (types[i], types[j], types[k])
        for (i, j, k), t in zip(idx, theta0) if np.isnan(t)
    )
    check_angles("GRO_BEFORE", gro_before, xyz_b, idx, theta0, mols, resnames, names, False)
    check_angles("GRO_AFTER", gro_after, xyz_a, idx, theta0, mols, resnames, names, True)
    if unparam:
        print("  unparameterised type triples (count):")
        for trip, cnt in unparam.most_common(MAX_EXAMPLES):
            print(f"    {'-'.join(trip)} {cnt}")

    # ---- Check B: ring planarity -------------------------------------
    print("\n== Check B: pyridinium ring planarity (NZ C10 C11 C12 C13 C14)")
    rows = []
    incomplete = 0
    for offset, molname in instances:
        incomplete += cache[molname]["incomplete"]
        for ring in cache[molname]["rings"]:
            rows.append([offset + i for i in ring])
    rings = np.array(rows, dtype=int).reshape(-1, 6)
    print(f"  incomplete ring residues skipped: {incomplete}")
    for label, gro, xyz in (("GRO_BEFORE", gro_before, xyz_b), ("GRO_AFTER", gro_after, xyz_a)):
        if not len(rings):
            print(f"{label}: {gro} rings=0")
            continue
        oop = ring_oop(xyz, rings)
        print(f"{label}: {gro} rings={len(oop)} mean={oop.mean():.3f} A max={oop.max():.3f} A "
              f"|>0.1 A: {(oop > 0.1).sum()} >0.3 A: {(oop > 0.3).sum()}")

    # ---- Check C: chirality flips ------------------------------------
    print("\n== Check C: chirality (signed volume) flips")
    adj = defaultdict(list)
    interest = {i for i in range(n_top) if names[i] == "CA" or resnames[i] in CHIRAL_RES}
    for offset, molname in instances:
        for a, b in cache[molname]["bonds"]:
            ga, gb = offset + a - 1, offset + b - 1
            if ga in interest:
                adj[ga].append(gb)
            if gb in interest:
                adj[gb].append(ga)

    centres, nbrs, is_ca = [], [], []
    for c in sorted(interest):
        if len(adj[c]) == 4:
            if sum(names[x].startswith("H") for x in adj[c]) >= 2:
                continue
            n1, n2, n3, _ = sorted(adj[c])
            centres.append(c)
            nbrs.append((n1, n2, n3))
            is_ca.append(names[c] == "CA")
    centres = np.array(centres, dtype=int)
    nbrs = np.array(nbrs, dtype=int).reshape(-1, 3)
    is_ca = np.array(is_ca, dtype=bool)
    print(f"  centres (4 bonded neighbours): CA={is_ca.sum()} "
          f"crosslink-residue other={(~is_ca).sum()}")

    v_b = signed_volumes(xyz_b, centres, nbrs)
    v_a = signed_volumes(xyz_a, centres, nbrs)
    flip = np.sign(v_b) != np.sign(v_a)
    print(f"  flipped GRO_BEFORE->GRO_AFTER: CA={(flip & is_ca).sum()} "
          f"crosslink-residue other={(flip & ~is_ca).sum()}")
    for c in centres[flip][:MAX_EXAMPLES]:
        print(f"    {mols[c]} {resnames[c]} {names[c]} {c}")

    # CA handedness in GRO_AFTER, non-GLY, neighbour set (N, C, CB) in that order
    ca_ng, missing = [], 0
    for c in sorted(interest):
        if names[c] != "CA" or resnames[c] == "GLY":
            continue
        byname = {names[x]: x for x in adj[c]}
        if all(k in byname for k in ("N", "C", "CB")):
            ca_ng.append((c, byname["N"], byname["C"], byname["CB"]))
        else:
            missing += 1
    if ca_ng:
        ca_ng = np.array(ca_ng, dtype=int)
        v_ca = signed_volumes(xyz_a, ca_ng[:, 0], ca_ng[:, 1:])
        print(f"  GRO_AFTER CA non-GLY (N,C,CB) signed volume: positive={(v_ca > 0).sum()} "
              f"negative={(v_ca < 0).sum()} missing_neighbour={missing}")


if __name__ == "__main__":
    main()
