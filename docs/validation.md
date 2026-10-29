# Validating a Topology

The scripts in [`scripts/validation/`](../scripts/validation/) check that a
built system has the expected crosslinks and runs in GROMACS. They need the
`colbuilder` environment (Python with numpy) and GROMACS (`gmx`).

## Build and check in one step

Put a `config.yaml` in an empty directory and run:

```bash
GMXRC=/path/to/gromacs/bin/GMXRC bash scripts/validation/run_topology_check.sh my_test/
```

This builds the system, then runs a short vacuum energy minimisation (EM),
solvates, adds ions and runs a short solvated EM. Afterwards it runs the
crosslink checks below. Each stage writes one line to `my_test/STATUS`, and the
last line is `DONE OK` or `DONE FAIL`:

```
BUILD OK 287s top=homo_sapiens_topology_files/collagen_fibril_homo_sapiens.top
EDITCONF OK 16.61700 18.58700 63.75302
GROMPP_VAC OK 1 warnings, 2 notes
EM_VAC OK Epot=-8.8871062e+05 Fmax=7.6162168e+03
SOLVATE OK 1932525 atoms (spc216.gro)
GENION OK Will try to add 1779 NA ions and 2061 CL ions.
GROMPP_SOLV OK 0 warnings, 0 notes
EM_SOLV OK Epot=-2.6199832e+07 Fmax=2.0819131e+04
UNBONDED OK close_pairs=20 bonded=20 unbonded=0 unbonded_same_molecule=0 unbonded_cross_molecule=0
XLINK_TERMS OK
XLINK_GEOMETRY OK
DONE OK 379s
```

Optional environment variables: `GMXRC` (a GROMACS environment script to
source), `GMX_THREADS` (OpenMP threads, default 8) and `VALIDATION_LOCK` (a
file to lock around the GROMACS steps when several checks run at once).

## Individual checks

Run these on any existing topology directory (the one holding
`collagen_fibril_*.top`, `.gro` and `col_*.itp`):

| Script | Usage | What to expect |
|--------|-------|----------------|
| `unbonded_xlinks.py` | `TOP GRO ITP_DIR` | Pairs of crosslink atoms on different molecules closer than 5 Å that are not bonded. Expect `unbonded=0`; the script exits 1 otherwise. |
| `xlink_terms.py` | `ITP_DIR [ITP_DIR2]` | The angles, dihedrals and 1-4 pairs that cross each crosslink bond. With two directories (for example a `legacy` and a `specbond` build of the same config), the `DIFF` section should say `no differences`. |
| `xlink_geometry.py` | `TOP ITP_DIR GRO [GRO ...]` | Crosslink bond lengths next to the Amber reference length. Bonds come out of the build stretched (typically 1.8-2.5 Å) and should relax to the reference (about 1.5 Å) after the vacuum EM. |
| `xlink_strain.py` | `TOP ITP_DIR FF_DIR GRO_BEFORE GRO_AFTER` | Angle strain around crosslinks, pyridinium ring planarity, and chirality flips between two structures. It also counts the handedness of every non-GLY CA, so D-residues show up as a small minority group. |

## Known issues

- **D-residues in PYD builds.** PYD fibrils can contain a D-configured LY2
  at 5.B (and sometimes LYX at 104.C). It comes from the sequence (MODELLER)
  step, not from topology generation. It shows up in the CA handedness
  count from `xlink_strain.py`.
- **Fibril clipping.** An experimental change to how the matrix set and crystal
  builder clip molecules at the fibril boundary is not included. Geometry is
  unchanged from upstream.
- **Splitting the final fibril (`specbond` mode).** Molecules are split where
  the chain ID changes; `TER` and `MODEL` records are not used. ColBuilder's own
  fibrils always change chain ID between chains. A hand-edited PDB in which two
  consecutive chains share an ID would be read as one chain.
- **Speed (`specbond` mode).** Crosslink bond verification compares each
  expected bond against every atom of the same name. This is fine for fibrils
  up to about 10^6 atoms, but will slow down for much larger components.
