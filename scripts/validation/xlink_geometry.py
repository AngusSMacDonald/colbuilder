"""Report crosslink bond lengths in one or more .gro files of the same system.

usage: python xlink_geometry.py TOP_FILE ITP_DIR GRO_FILE [GRO_FILE ...]

Crosslink bonds are the divalent CE-NZ, CE-NZ1 and CE1-CE1 bonds and the
trivalent C12-CB, C13-CG and C12-CG bonds in the topology, between two
different residues. Distances are printed in Angstrom next to the Amber
equilibrium length, so a structure built with stretched crosslinks shows up
before and after minimisation.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unbonded_xlinks import parse_itp, parse_top_molecules  # noqa: E402

# (atom name pair) -> (label, equilibrium length in A from ffbonded.itp)
XLINK_BONDS = {
    frozenset(("CE", "NZ")): ("CE-NZ (divalent)", 1.471),   # CT-N3
    frozenset(("CE", "NZ1")): ("CE-NZ1 (MOLD)", 1.471),
    frozenset(("C12", "CB")): ("C12-CB (trivalent)", 1.510),  # CA-CT
    frozenset(("C13", "CG")): ("C13-CG (trivalent)", 1.510),  # CA-CT
    frozenset(("C12", "CG")): ("C12-CG (PYL/DPL)", 1.510),
    frozenset(("CE1",)): ("CE1-CE1 (DTY-DTF)", 1.490),  # both ends named CE1
}


def crosslink_bonds(top_file, itp_dir):
    found, offset = [], 0
    for molname, count in parse_top_molecules(top_file):
        atoms, bonds = parse_itp(os.path.join(itp_dir, molname + ".itp"))
        for _ in range(count):
            for a, b in bonds:
                key = frozenset((atoms[a - 1][1], atoms[b - 1][1]))
                if key in XLINK_BONDS and atoms[a - 1][0] != atoms[b - 1][0]:
                    found.append((XLINK_BONDS[key], offset + a - 1, offset + b - 1))
            offset += len(atoms)
    return found


def read_coords(gro_file, n):
    with open(gro_file) as fh:
        fh.readline()
        fh.readline()
        xyz = [fh.readline() for _ in range(n)]
    return np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])] for l in xyz]) * 10.0


def main():
    top_file, itp_dir, gros = sys.argv[1], sys.argv[2], sys.argv[3:]
    bonds = crosslink_bonds(top_file, itp_dir)
    if not bonds:
        print("no crosslink bonds in topology")
        return
    n = max(max(a, b) for _, a, b in bonds) + 1
    for gro in gros:
        xyz = read_coords(gro, n)
        print(f"== {gro}")
        for label, b0 in sorted({lab for lab, _, _ in bonds}):
            d = np.array([np.linalg.norm(xyz[a] - xyz[b]) for (lab, _), a, b in bonds if lab == label])
            print(f"  {label:16} n={len(d):4}  b0={b0:.2f}  min={d.min():.2f}  "
                  f"mean={d.mean():.2f}  max={d.max():.2f}  >1.7A: {(d > 1.7).sum()}")


if __name__ == "__main__":
    main()
