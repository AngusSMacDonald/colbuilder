#!/usr/bin/env python3
"""Count geometrically close crosslink atom pairs that are not covalently bonded.

Usage: python3 unbonded_xlinks.py TOP_FILE GRO_FILE ITP_DIR

Exits with status 1 if any close crosslink pair is unbonded.
"""
import os
import re
import sys
from collections import defaultdict

import numpy as np

CUTOFF_NM = 0.5
MAX_EXAMPLES = 10

# Compatible crosslink marker pairs: (resname, atomname) <-> (resname, atomname)
XLINK_PAIRS = [
    # trivalent
    (("LYX", "C12"), ("LY2", "CB")),
    (("LYX", "C13"), ("LY3", "CG")),
    (("LXX", "C12"), ("LX2", "CB")),
    (("LXX", "C13"), ("LX3", "CG")),
    (("LXY", "C12"), ("L2Y", "CG")),
    (("LXY", "C13"), ("L3Y", "CG")),
    (("LYY", "C12"), ("L2X", "CG")),
    (("LYY", "C13"), ("L3X", "CG")),
    # divalent
    (("L4Y", "CE"), ("L5Y", "NZ")),
    (("L4X", "CE"), ("L5X", "NZ")),
    (("LY4", "CE"), ("LY5", "NZ")),
    (("LX4", "CE"), ("LX5", "NZ")),
    (("LGX", "CE"), ("AGS", "NZ")),
    (("LPS", "CE"), ("APD", "NZ")),
    (("LZD", "CE"), ("LZS", "NZ1")),
    (("DTY", "CE1"), ("DTF", "CE1")),
]

SECTION_RE = re.compile(r"^\s*\[\s*(\S+?)\s*\]")


def strip_comment(line):
    return line.split(";", 1)[0].strip()


def parse_top_molecules(top_file):
    """Return list of (molname, count) from the [ molecules ] section."""
    molecules = []
    in_section = False
    with open(top_file) as fh:
        for raw in fh:
            m = SECTION_RE.match(raw)
            if m:
                in_section = (m.group(1) == "molecules")
                continue
            if not in_section:
                continue
            line = strip_comment(raw)
            if not line:
                continue
            parts = line.split()
            molecules.append((parts[0], int(parts[1])))
    return molecules


def parse_itp(itp_file):
    """Return (atoms, bonds) for one ITP file.

    atoms: list of (resname, atomname), index = local nr - 1
    bonds: list of (local_a, local_b), 1-based local atom numbers
    Stops reading at the [ pairs ] section.
    """
    atoms_by_nr = {}
    bonds = []
    section = None
    with open(itp_file) as fh:
        for raw in fh:
            if raw[:1] == "[":
                m = SECTION_RE.match(raw)
                if m:
                    section = m.group(1)
                    if section == "pairs":
                        break
                    continue
            if raw[:1] == ";":
                continue
            line = strip_comment(raw)
            if not line:
                continue
            if section == "atoms":
                parts = line.split()
                # nr type resnr residue atom cgnr charge mass
                atoms_by_nr[int(parts[0])] = (parts[3], parts[4])
            elif section == "bonds":
                parts = line.split()
                bonds.append((int(parts[0]), int(parts[1])))
    atoms = [atoms_by_nr[nr] for nr in sorted(atoms_by_nr)]
    return atoms, bonds


def main():
    if len(sys.argv) != 4:
        print("usage: python3 unbonded_xlinks.py TOP_FILE GRO_FILE ITP_DIR", file=sys.stderr)
        sys.exit(2)
    top_file, gro_file, itp_dir = sys.argv[1], sys.argv[2], sys.argv[3]

    # ---- Topology: build global atom list and bond set -----------------
    molecules = parse_top_molecules(top_file)
    itp_cache = {}
    global_res = []
    global_name = []
    global_mol = []      # molecule name per global atom
    global_inst = []     # molecule instance id per global atom
    bond_set = set()
    skipped = []
    inst_id = 0

    for molname, count in molecules:
        itp_file = os.path.join(itp_dir, molname + ".itp")
        if not os.path.isfile(itp_file):
            skipped.append(molname)
            continue
        if molname not in itp_cache:
            itp_cache[molname] = parse_itp(itp_file)
        atoms, bonds = itp_cache[molname]
        n_local = len(atoms)
        for _ in range(count):
            offset = len(global_res)
            for res, name in atoms:
                global_res.append(res)
                global_name.append(name)
                global_mol.append(molname)
                global_inst.append(inst_id)
            for a, b in bonds:
                bond_set.add(frozenset((offset + a - 1, offset + b - 1)))
            inst_id += 1
            assert len(global_res) == offset + n_local

    if skipped:
        print("WARNING: no ITP file for molecules (ignored): " + ", ".join(sorted(set(skipped))))

    n_top = len(global_res)

    # ---- GRO: coordinates, line order as global index -----------------
    with open(gro_file) as fh:
        fh.readline()  # title line
        n_gro = int(fh.readline().strip())
        gro_res = []
        gro_name = []
        coords = np.empty((n_gro, 3), dtype=float)
        for i in range(n_gro):
            line = fh.readline()
            gro_res.append(line[5:10].strip())
            gro_name.append(line[10:15].strip())
            coords[i, 0] = float(line[20:28])
            coords[i, 1] = float(line[28:36])
            coords[i, 2] = float(line[36:44])

    # ---- Sanity checks ------------------------------------------------
    if n_gro != n_top:
        print(f"SANITY FAIL: gro atoms={n_gro} != itp atoms summed over molecules={n_top}")
        sys.exit(1)
    mismatches = []
    for i in range(n_gro):
        if gro_res[i] != global_res[i] or gro_name[i] != global_name[i]:
            mismatches.append(
                f"#{i + 1}: gro {gro_res[i]}:{gro_name[i]} vs itp "
                f"{global_res[i]}:{global_name[i]} ({global_mol[i]})"
            )
    if mismatches:
        print(f"SANITY FAIL: {len(mismatches)} gro/itp atom mismatches; first 5:")
        for msg in mismatches[:5]:
            print("  " + msg)
        sys.exit(1)

    # ---- Marker atom indices per type --------------------------------
    by_type = defaultdict(list)
    for i, key in enumerate(zip(global_res, global_name)):
        by_type[key].append(i)

    # ---- Close pairs per compatible type pair -------------------------
    results = []  # (typeA, typeB, idxA, idxB, dist_nm, bonded)
    per_type = []  # (labelA, labelB, n_bonded, n_unbonded)
    for ta, tb in XLINK_PAIRS:
        idx_a = np.array(by_type.get(ta, []), dtype=int)
        idx_b = np.array(by_type.get(tb, []), dtype=int)
        n_bonded = 0
        n_unbonded = 0
        if idx_a.size and idx_b.size:
            pa = coords[idx_a]
            pb = coords[idx_b]
            diff = pa[:, None, :] - pb[None, :, :]
            dist = np.sqrt((diff * diff).sum(axis=2))
            ii, jj = np.nonzero(dist < CUTOFF_NM)
            for i, j in zip(ii, jj):
                ga = int(idx_a[i])
                gb = int(idx_b[j])
                bonded = frozenset((ga, gb)) in bond_set
                results.append((ta, tb, ga, gb, float(dist[i, j]), bonded))
                if bonded:
                    n_bonded += 1
                else:
                    n_unbonded += 1
        per_type.append((ta, tb, n_bonded, n_unbonded))

    results.sort(key=lambda r: (r[2], r[3]))
    unbonded = [r for r in results if not r[5]]
    n_close = len(results)
    n_bonded_all = n_close - len(unbonded)
    n_same = 0
    n_cross = 0
    for _, _, ga, gb, _, _ in unbonded:
        if global_inst[ga] == global_inst[gb]:
            n_same += 1
        else:
            n_cross += 1

    label = os.path.dirname(os.path.abspath(top_file))
    print(
        f"{label} close_pairs={n_close} bonded={n_bonded_all} "
        f"unbonded={len(unbonded)} unbonded_same_molecule={n_same} "
        f"unbonded_cross_molecule={n_cross}"
    )

    for ta, tb, ga, gb, d, _ in unbonded[:MAX_EXAMPLES]:
        print(
            f"UNBONDED {global_mol[ga]} res{ta[0]}:{ta[1]} #{ga + 1} <-> "
            f"{global_mol[gb]} {tb[0]}:{tb[1]} #{gb + 1} {d * 10:.2f} A"
        )

    print("PER_TYPE breakdown:")
    for ta, tb, nb, nu in per_type:
        print(f"  {ta[0]}:{ta[1]} <-> {tb[0]}:{tb[1]} bonded={nb} unbonded={nu}")

    if unbonded:
        sys.exit(1)


if __name__ == "__main__":
    main()
