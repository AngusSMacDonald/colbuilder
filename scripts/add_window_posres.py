#!/usr/bin/env python3
"""
Add window-based end position restraints to Martini collagen fibrils.

This helper reads a final ColBuilder Martini topology directory containing:
  - collagen_fibril_CG_*.pdb
  - col_<n>.itp

For each `col_<n>.itp`, it:
  1. maps the molecule's atoms onto the corresponding slice of the combined CG PDB,
  2. identifies candidate beads (`BB` by default) near one or both ends along an axis,
  3. rewrites the `[ position_restraints ]` section to contain only those atoms.

Selection can be done either by geometric end window (`--window-nm`) or by
an exact candidate count per molecule end (`--terminal-count`).
"""

from __future__ import annotations

import argparse
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence


ITP_PATTERN = re.compile(r"^col_(\d+)\.itp$")
SECTION_HEADER_PATTERN = re.compile(r"^\[\s*(.+?)\s*\]$")


@dataclass(frozen=True)
class AtomRecord:
    nr: int
    atom_name: str


@dataclass(frozen=True)
class CoordRecord:
    x: float
    y: float
    z: float

    def value(self, axis: str) -> float:
        return getattr(self, axis)


@dataclass(frozen=True)
class PdbAtomRecord:
    index: int
    serial: int
    atom_name: str
    chain_id: str
    coord: CoordRecord
    raw_line: str


@dataclass
class MoleculeSummary:
    itp_name: str
    total_atoms: int
    candidate_atoms: int
    restrained_atoms: int
    min_coord: float
    max_coord: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Restrict Martini position restraints to fibril-end windows."
    )
    parser.add_argument(
        "--topology-dir",
        type=Path,
        default=Path.cwd(),
        help="Directory containing collagen_fibril_CG_*.pdb and col_<n>.itp files.",
    )
    parser.add_argument(
        "--pdb",
        type=Path,
        help="Combined CG PDB file. If omitted, auto-detect collagen_fibril_CG_*.pdb in --topology-dir.",
    )
    parser.add_argument(
        "--axis",
        choices=("x", "y", "z"),
        default="z",
        help="Axis used to define low/high fibril ends.",
    )
    parser.add_argument(
        "--end",
        choices=("both", "low", "high"),
        default="both",
        help="Which fibril end to restrain.",
    )
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "--window-nm",
        type=float,
        help="End-window thickness in nm. The script converts this to Angstrom for the PDB coordinates.",
    )
    mode_group.add_argument(
        "--terminal-count",
        type=int,
        help="Select the exact number of candidate atoms at each requested end of each molecule.",
    )
    parser.add_argument(
        "--atom-name",
        default="BB",
        help="Atom/bead name to restrain. Default: BB.",
    )
    parser.add_argument(
        "--force",
        nargs="+",
        type=str,
        metavar="FC",
        default=("1000", "1000", "1000"),
        help=(
            "Position-restraint force constants. Provide either one value "
            "(for example POSRES_FC or 1000) to use for x/y/z, or three values "
            "for separate x/y/z entries."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report the selection without modifying any ITP files.",
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not copy original ITP files to a backup directory before rewriting.",
    )
    parser.add_argument(
        "--backup-dir",
        type=Path,
        default=Path("backup"),
        help="Backup directory written relative to --topology-dir unless absolute.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("end_posre_summary.tsv"),
        help="Summary file written relative to --topology-dir. Use '-' to skip writing a summary.",
    )
    parser.add_argument(
        "--restrained-pdb",
        type=Path,
        help="Optional output PDB containing only the restrained atoms for easy VMD loading.",
    )
    parser.add_argument(
        "--vmd-tcl",
        type=Path,
        help="Optional VMD Tcl script defining an exact restrained-atom selection on the full CG PDB.",
    )
    return parser.parse_args()


def parse_force_constants(values: Sequence[str]) -> List[str]:
    if len(values) == 1:
        return [values[0], values[0], values[0]]
    if len(values) == 3:
        return list(values)
    raise ValueError("--force expects either one value or three values")


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def autodetect_pdb(topology_dir: Path) -> Path:
    candidates = sorted(topology_dir.glob("collagen_fibril_CG_*.pdb"))
    if not candidates:
        raise FileNotFoundError(
            f"No collagen_fibril_CG_*.pdb file found in {topology_dir}"
        )
    if len(candidates) > 1:
        joined = ", ".join(item.name for item in candidates)
        raise ValueError(
            f"Multiple collagen_fibril_CG_*.pdb files found in {topology_dir}: {joined}. "
            "Pass --pdb explicitly."
        )
    return candidates[0].resolve()


def discover_itps(topology_dir: Path) -> List[Path]:
    matches = []
    for path in topology_dir.glob("col_*.itp"):
        match = ITP_PATTERN.match(path.name)
        if match:
            matches.append((int(match.group(1)), path.resolve()))
    if not matches:
        raise FileNotFoundError(f"No col_<n>.itp files found in {topology_dir}")
    return [path for _, path in sorted(matches)]


def parse_itp_atoms(itp_path: Path) -> List[AtomRecord]:
    atoms: List[AtomRecord] = []
    current_section = None

    with itp_path.open() as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped or stripped.startswith(";"):
                continue

            section_match = SECTION_HEADER_PATTERN.match(stripped)
            if section_match:
                current_section = section_match.group(1).lower()
                continue

            if current_section != "atoms" or stripped.startswith("#"):
                continue

            fields = stripped.split()
            if len(fields) < 5:
                continue

            try:
                atoms.append(AtomRecord(nr=int(fields[0]), atom_name=fields[4]))
            except ValueError as exc:
                raise ValueError(
                    f"Failed parsing atom line in {itp_path}: {raw_line.rstrip()}"
                ) from exc

    if not atoms:
        raise ValueError(f"No [ atoms ] entries found in {itp_path}")
    return atoms


def parse_pdb_atoms(pdb_path: Path) -> List[PdbAtomRecord]:
    atoms: List[PdbAtomRecord] = []

    with pdb_path.open() as handle:
        for raw_line in handle:
            if not raw_line.startswith(("ATOM", "HETATM")):
                continue

            try:
                serial = int(raw_line[6:11].strip())
                atom_name = raw_line[12:16].strip()
                chain_id = raw_line[21:22].strip() or "_"
                coord = CoordRecord(
                    x=float(raw_line[30:38].strip()),
                    y=float(raw_line[38:46].strip()),
                    z=float(raw_line[46:54].strip()),
                )
            except ValueError as exc:
                raise ValueError(
                    f"Failed parsing PDB atom line in {pdb_path}: {raw_line.rstrip()}"
                ) from exc

            atoms.append(
                PdbAtomRecord(
                    index=len(atoms),
                    serial=serial,
                    atom_name=atom_name,
                    chain_id=chain_id,
                    coord=coord,
                    raw_line=raw_line.rstrip("\n"),
                )
            )

    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {pdb_path}")
    return atoms


def chunk_records(
    records: Sequence[PdbAtomRecord], atom_counts: Iterable[int]
) -> List[List[PdbAtomRecord]]:
    chunks: List[List[PdbAtomRecord]] = []
    start = 0
    for count in atom_counts:
        end = start + count
        chunk = list(records[start:end])
        if len(chunk) != count:
            raise ValueError(
                "Combined PDB atom count does not match the sum of atoms in the ITP files"
            )
        chunks.append(chunk)
        start = end

    if start != len(records):
        raise ValueError(
            "Combined PDB contains extra atoms beyond the sum of the discovered ITP files"
        )
    return chunks


def select_end_atoms(
    atoms: Sequence[AtomRecord],
    pdb_atoms: Sequence[PdbAtomRecord],
    axis: str,
    atom_name: str,
    end_mode: str,
    window_nm: float | None,
    terminal_count: int | None,
) -> tuple[MoleculeSummary, List[int], List[PdbAtomRecord]]:
    if len(atoms) != len(pdb_atoms):
        raise ValueError("Atom and coordinate counts do not match for molecule selection")

    candidates = []
    for atom, pdb_atom in zip(atoms, pdb_atoms):
        if atom.atom_name != atom_name:
            continue
        candidates.append((pdb_atom.coord.value(axis), atom.nr, pdb_atom))

    if not candidates:
        raise ValueError(f"No atoms named '{atom_name}' were found in the molecule")

    coords = [value for value, _, _ in candidates]
    min_coord = min(coords)
    max_coord = max(coords)

    selected_atom_ids: set[int] = set()

    if window_nm is not None:
        window_a = window_nm * 10.0
        low_cutoff = min_coord + window_a
        high_cutoff = max_coord - window_a

        if end_mode in {"both", "low"}:
            for coord_value, atom_nr, _ in candidates:
                if coord_value <= low_cutoff:
                    selected_atom_ids.add(atom_nr)

        if end_mode in {"both", "high"}:
            for coord_value, atom_nr, _ in candidates:
                if coord_value >= high_cutoff:
                    selected_atom_ids.add(atom_nr)
    else:
        if terminal_count is None:
            raise ValueError("Either --window-nm or --terminal-count must be provided")

        ordered = sorted(candidates, key=lambda item: (item[0], item[1]))
        n_select = min(terminal_count, len(ordered))

        if end_mode in {"both", "low"}:
            for _, atom_nr, _ in ordered[:n_select]:
                selected_atom_ids.add(atom_nr)

        if end_mode in {"both", "high"}:
            for _, atom_nr, _ in ordered[-n_select:]:
                selected_atom_ids.add(atom_nr)

    restrained_atoms: List[int] = []
    restrained_pdb_atoms: List[PdbAtomRecord] = []
    for atom, pdb_atom in zip(atoms, pdb_atoms):
        if atom.nr in selected_atom_ids:
            restrained_atoms.append(atom.nr)
            restrained_pdb_atoms.append(pdb_atom)

    summary = MoleculeSummary(
        itp_name="",
        total_atoms=len(atoms),
        candidate_atoms=len(candidates),
        restrained_atoms=len(restrained_atoms),
        min_coord=min_coord,
        max_coord=max_coord,
    )
    return summary, restrained_atoms, restrained_pdb_atoms


def rewrite_position_restraints(
    itp_path: Path,
    restrained_atoms: Sequence[int],
    force_constants: Sequence[str],
    backup_dir: Path | None,
) -> None:
    lines = itp_path.read_text().splitlines(keepends=True)

    if backup_dir is not None:
        backup_dir.mkdir(parents=True, exist_ok=True)
        backup_path = backup_dir / itp_path.name
        if not backup_path.exists():
            shutil.copy2(itp_path, backup_path)

    block = ["[ position_restraints ]\n", "#ifdef POSRES\n"]
    for atom_nr in restrained_atoms:
        block.append(
            f"{atom_nr} 1 {force_constants[0]} {force_constants[1]} {force_constants[2]}\n"
        )
    block.append("#endif\n")

    section_start = None
    section_end = None

    for index, line in enumerate(lines):
        stripped = line.strip()
        section_match = SECTION_HEADER_PATTERN.match(stripped)
        if section_start is None:
            if section_match and section_match.group(1).lower() == "position_restraints":
                section_start = index
        else:
            if section_match:
                section_end = index
                break

    if section_start is None:
        new_lines = list(lines)
        if new_lines and new_lines[-1].strip():
            new_lines.append("\n")
        new_lines.extend(block)
    else:
        if section_end is None:
            section_end = len(lines)
        new_lines = lines[:section_start] + block + lines[section_end:]

    itp_path.write_text("".join(new_lines))


def write_summary(summary_path: Path, summaries: Sequence[MoleculeSummary]) -> None:
    ensure_parent(summary_path)
    with summary_path.open("w") as handle:
        handle.write(
            "itp\ttotal_atoms\tcandidate_atoms\trestrained_atoms\tmin_coord_A\tmax_coord_A\n"
        )
        for item in summaries:
            handle.write(
                f"{item.itp_name}\t{item.total_atoms}\t{item.candidate_atoms}\t"
                f"{item.restrained_atoms}\t{item.min_coord:.3f}\t{item.max_coord:.3f}\n"
            )


def write_restrained_pdb(
    restrained_pdb_path: Path, restrained_pdb_atoms: Sequence[PdbAtomRecord]
) -> None:
    ensure_parent(restrained_pdb_path)
    with restrained_pdb_path.open("w") as handle:
        for atom in restrained_pdb_atoms:
            handle.write(f"{atom.raw_line}\n")
        handle.write("END\n")


def write_vmd_tcl(
    vmd_tcl_path: Path, indices: Sequence[int], expected_num_atoms: int
) -> None:
    ensure_parent(vmd_tcl_path)
    joined_indices = " ".join(str(index) for index in indices)
    max_index = max(indices) if indices else -1
    with vmd_tcl_path.open("w") as handle:
        handle.write("# Auto-generated by add_window_posres.py\n")
        handle.write(f"set posre_expected_numatoms {expected_num_atoms}\n")
        handle.write(f"set posre_expected_selected {len(indices)}\n")
        handle.write(f"set posre_max_index {max_index}\n")
        handle.write(f"set posre_end_indices {{{joined_indices}}}\n")
        handle.write('set posre_top_numatoms [molinfo top get numatoms]\n')
        handle.write('if {$posre_top_numatoms <= $posre_max_index} {\n')
        handle.write(
            '    puts "ERROR: top molecule has only $posre_top_numatoms atoms, but the '
            'selection needs index $posre_max_index. Load the matching system or regenerate '
            'the selection."\n'
        )
        handle.write('}\n')
        handle.write(
            'if {$posre_top_numatoms > $posre_expected_numatoms} {\n'
            '    puts "NOTE: top molecule has extra atoms beyond the collagen-only system '
            '(likely solvent/ions). Selection assumes the collagen atom ordering is unchanged."\n'
            '}\n'
        )
        handle.write('set posre_end_atoms [atomselect top "index $posre_end_indices"]\n')
        handle.write('if {[$posre_end_atoms num] != $posre_expected_selected} {\n')
        handle.write(
            '    puts "WARNING: selected [$posre_end_atoms num] atoms, expected '
            '$posre_expected_selected. Check atom ordering and the loaded molecule."\n'
        )
        handle.write('}\n')
        handle.write('mol addrep top\n')
        handle.write('set posre_rep [expr {[molinfo top get numreps] - 1}]\n')
        handle.write('mol modselect $posre_rep top "index $posre_end_indices"\n')
        handle.write('mol modstyle $posre_rep top VDW 1.0 12.0\n')
        handle.write('mol modcolor $posre_rep top ColorID 1\n')


def main() -> None:
    args = parse_args()
    if args.window_nm is not None and args.window_nm <= 0:
        raise ValueError("--window-nm must be positive")
    if args.terminal_count is not None and args.terminal_count <= 0:
        raise ValueError("--terminal-count must be positive")

    force_constants = parse_force_constants(args.force)

    topology_dir = args.topology_dir.resolve()
    pdb_path = args.pdb.resolve() if args.pdb else autodetect_pdb(topology_dir)
    itp_paths = discover_itps(topology_dir)

    backup_dir = None
    if not args.no_backup:
        backup_dir = args.backup_dir
        if not backup_dir.is_absolute():
            backup_dir = topology_dir / backup_dir

    summary_path = None
    if args.summary != Path("-"):
        summary_path = args.summary
        if not summary_path.is_absolute():
            summary_path = topology_dir / summary_path

    restrained_pdb_path = None
    if args.restrained_pdb:
        restrained_pdb_path = args.restrained_pdb
        if not restrained_pdb_path.is_absolute():
            restrained_pdb_path = topology_dir / restrained_pdb_path

    vmd_tcl_path = None
    if args.vmd_tcl:
        vmd_tcl_path = args.vmd_tcl
        if not vmd_tcl_path.is_absolute():
            vmd_tcl_path = topology_dir / vmd_tcl_path

    pdb_atoms = parse_pdb_atoms(pdb_path)
    parsed_itps = [(path, parse_itp_atoms(path)) for path in itp_paths]
    pdb_chunks = chunk_records(
        pdb_atoms,
        atom_counts=(len(atoms) for _, atoms in parsed_itps),
    )

    summaries: List[MoleculeSummary] = []
    all_restrained_pdb_atoms: List[PdbAtomRecord] = []

    for (itp_path, atoms), pdb_chunk in zip(parsed_itps, pdb_chunks):
        summary, restrained_atoms, restrained_pdb_atoms = select_end_atoms(
            atoms=atoms,
            pdb_atoms=pdb_chunk,
            axis=args.axis,
            atom_name=args.atom_name,
            end_mode=args.end,
            window_nm=args.window_nm,
            terminal_count=args.terminal_count,
        )
        summary.itp_name = itp_path.name
        summaries.append(summary)
        all_restrained_pdb_atoms.extend(restrained_pdb_atoms)

        if not args.dry_run:
            rewrite_position_restraints(
                itp_path=itp_path,
                restrained_atoms=restrained_atoms,
                force_constants=force_constants,
                backup_dir=backup_dir,
            )

    if summary_path is not None:
        write_summary(summary_path, summaries)

    if restrained_pdb_path is not None:
        write_restrained_pdb(restrained_pdb_path, all_restrained_pdb_atoms)

    if vmd_tcl_path is not None:
        write_vmd_tcl(
            vmd_tcl_path,
            indices=[atom.index for atom in all_restrained_pdb_atoms],
            expected_num_atoms=len(pdb_atoms),
        )

    action = "Would rewrite" if args.dry_run else "Rewrote"
    selection_mode = (
        f"window={args.window_nm} nm"
        if args.window_nm is not None
        else f"terminal-count={args.terminal_count}"
    )
    print(
        f"{action} {len(summaries)} ITP files in {topology_dir} using "
        f"axis={args.axis}, {selection_mode}, end={args.end}, atom={args.atom_name}."
    )
    for item in summaries:
        print(
            f"{item.itp_name}: restrained {item.restrained_atoms}/{item.candidate_atoms} "
            f"candidate atoms between {item.min_coord:.3f} A and {item.max_coord:.3f} A"
        )


if __name__ == "__main__":
    main()
