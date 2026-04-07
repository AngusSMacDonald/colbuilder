These are helper scripts to build topology for mixed fibrils. You need to have run pdb (triple-helix) generation for all crosslink types required first, and then the mixing script with force field specified but topology generation off.

i.e. sample workflow in colbuilder/colbuilder-dev environment
1. triple_helix_HLKNL.yaml (sequence generate w/ crosslink HLKNL at desired position, no geometry or topology generator, debug: true)
2. triple_helix_PYD.yaml (repeat for crosslink PYD, debug: true)
3. repeat _ad infinitum_ for however many crosslink types you want
3. mix_geometry.yaml 
 - mix_bool with ratio of two/more input pdb from above
 - set contact_distance and fibril_length as desired
 - set force_field: "martini3"
 - sequence/geometry/topology_generator: false
4. For Martini3:
 - `python scripts/topology_from_mix.py mix_geometry.yaml`
5. For atomistic AMBER99:
 - set `force_field: "amber99"` in the same mix config
 - `python scripts/topology_from_mix_atomistic.py mix_geometry.yaml`

Martini end restraints:
- Window or whole-fibril end slab mode: `python scripts/add_window_posres.py --topology-dir /path/to/homo_sapiens_martini3_topology_files --axis z --window-nm 2.0 --end both --dry-run`
- Exact chain-aware termini mode: `python scripts/add_termini_posres.py --topology-dir /path/to/homo_sapiens_martini3_topology_files --axis z --end high --per-chain-count 1 --dry-run`
- `add_window_posres.py` and `add_termini_posres.py` are standalone scripts. `filter_martini_end_posre.py` is now only a backward-compatible alias to the window-based script.
- Remove `--dry-run` to rewrite the existing `[ position_restraints ]` blocks in each `col_<n>.itp`.
- Original `col_<n>.itp` files are copied to `backup/` inside the topology directory by default.
- `--force` accepts either one value, e.g. `--force POSRES_FC`, or three values for separate x/y/z entries.
- `add_window_posres.py --terminal-count 1 --end high` selects exactly one BB at the high end; `--terminal-count 1 --end both` selects one BB at each end.
- `add_termini_posres.py --per-chain-count 1 --end high` selects one terminal BB per collagen chain at the high end.
- `add_termini_posres.py --index-file pull_groups.ndx` writes GROMACS pull groups:
  - by default only `TERM_LOW`, `TERM_HIGH`, `TERM_BOTH`
  - add `--index-detail molecule` for `TERM_COL0_LOW`, `TERM_COL0_HIGH`, etc.
  - add `--index-detail chain` for `TERM_COL0_A_LOW`, `TERM_COL0_B_HIGH`, etc.
  - add `--index-detail all` for every level
- `add_termini_posres.py --pull-mdp-snippet pull_groups.mdpfrag` writes a matching pull-code `.mdp` block. Most useful with `--index-detail molecule` or `--index-detail all`.
- Add `--merge-index` to preserve existing groups in an input `.ndx` and replace only any exported termini-group names that already exist.
- For VMD, you can also write `--restrained-pdb end_posre_atoms.pdb` and/or `--vmd-tcl end_posre_atoms.tcl`.
- These helpers are intended for tendon-like boundary conditions where only BB beads near the fibril ends should remain restrained.
