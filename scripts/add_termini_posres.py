#!/usr/bin/env python3
"""
Add exact chain-aware terminal position restraints to Martini collagen fibrils.

This helper reads a final ColBuilder Martini topology directory containing:
  - collagen_fibril_CG_*.pdb
  - col_<n>.itp

For each `col_<n>.itp`, it identifies backbone beads by chain and selects the
requested number of terminal BB beads at the low and/or high end of each
collagen chain.
"""

from __future__ import annotations

import argparse
import re
import shutil
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


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
class MoleculeTerminiSummary:
    itp_name: str
    total_atoms: int
    candidate_atoms: int
    chain_count: int
    restrained_atoms: int
    min_coord: float
    max_coord: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Restrict Martini position restraints to exact chain termini."
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
        help="Axis used to define the termini ordering.",
    )
    parser.add_argument(
        "--end",
        choices=("both", "low", "high"),
        default="high",
        help="Which terminus to restrain for each chain.",
    )
    parser.add_argument(
        "--per-chain-count",
        type=int,
        default=1,
        help="Number of candidate atoms to select at each requested end of each chain.",
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
        default=Path("termini_posre_summary.tsv"),
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
    parser.add_argument(
        "--index-file",
        type=Path,
        help="Optional GROMACS .ndx file containing termini pull groups.",
    )
    parser.add_argument(
        "--merge-index",
        action="store_true",
        help="Merge exported pull groups into an existing .ndx file instead of overwriting it.",
    )
    parser.add_argument(
        "--index-prefix",
        default="TERM",
        help="Prefix used for exported .ndx group names. Default: TERM.",
    )
    parser.add_argument(
        "--index-detail",
        choices=("aggregate", "molecule", "chain", "all"),
        default="aggregate",
        help=(
            "How many pull-group levels to export to the .ndx file. "
            "'aggregate' writes only PREFIX_LOW/HIGH/BOTH. "
            "'molecule' adds PREFIX_COLn_LOW/HIGH. "
            "'chain' adds PREFIX_COLn_A_LOW/HIGH-style groups. "
            "'all' writes every level."
        ),
    )
    parser.add_argument(
        "--pull-mdp-snippet",
        type=Path,
        help="Optional output file containing a matching GROMACS pull-code .mdp snippet.",
    )
    parser.add_argument(
        "--pull-type",
        choices=("umbrella", "constant-force"),
        default="umbrella",
        help="Pull-coordinate type written to --pull-mdp-snippet. Default: umbrella.",
    )
    parser.add_argument(
        "--pull-k",
        default="1000",
        help=(
            "Value written as pull-coord<n>-k in --pull-mdp-snippet. "
            "For constant-force, pass the force in kJ mol^-1 nm^-1, typically negative for extension."
        ),
    )
    parser.add_argument(
        "--pull-rate",
        default="0.0",
        help="Value written as pull-coord<n>-rate for umbrella snippets. Default: 0.0.",
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
        handle.write("# Auto-generated by add_termini_posres.py\n")
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


def chain_id_from_pdb_atom(pdb_atom: PdbAtomRecord) -> str:
    return pdb_atom.chain_id or "_"


def select_chain_termini(
    atoms: Sequence[AtomRecord],
    pdb_atoms: Sequence[PdbAtomRecord],
    axis: str,
    atom_name: str,
    end_mode: str,
    per_chain_count: int,
) -> tuple[
    MoleculeTerminiSummary,
    List[int],
    List[PdbAtomRecord],
    Dict[str, Dict[str, List[PdbAtomRecord]]],
]:
    if len(atoms) != len(pdb_atoms):
        raise ValueError("Atom and coordinate counts do not match for molecule selection")

    candidates_by_chain: Dict[str, List[tuple[float, int, PdbAtomRecord]]] = defaultdict(list)

    for atom, pdb_atom in zip(atoms, pdb_atoms):
        if atom.atom_name != atom_name:
            continue
        chain_id = chain_id_from_pdb_atom(pdb_atom)
        coord_value = pdb_atom.coord.value(axis)
        candidates_by_chain[chain_id].append((coord_value, atom.nr, pdb_atom))

    if not candidates_by_chain:
        raise ValueError(f"No atoms named '{atom_name}' were found in the molecule")

    all_coords = [value for items in candidates_by_chain.values() for value, _, _ in items]
    min_coord = min(all_coords)
    max_coord = max(all_coords)

    selected_atom_ids: set[int] = set()
    selected_by_end: Dict[str, Dict[str, List[PdbAtomRecord]]] = {"low": {}, "high": {}}
    for chain_id in sorted(candidates_by_chain):
        items = sorted(candidates_by_chain[chain_id], key=lambda item: (item[0], item[1]))
        n_select = min(per_chain_count, len(items))

        if end_mode in {"both", "low"}:
            low_atoms = [pdb_atom for _, _, pdb_atom in items[:n_select]]
            selected_by_end["low"][chain_id] = low_atoms
            for _, atom_nr, _ in items[:n_select]:
                selected_atom_ids.add(atom_nr)

        if end_mode in {"both", "high"}:
            high_atoms = [pdb_atom for _, _, pdb_atom in items[-n_select:]]
            selected_by_end["high"][chain_id] = high_atoms
            for _, atom_nr, _ in items[-n_select:]:
                selected_atom_ids.add(atom_nr)

    selected_atom_nrs: List[int] = []
    selected_pdb_atoms: List[PdbAtomRecord] = []
    for atom, pdb_atom in zip(atoms, pdb_atoms):
        if atom.nr in selected_atom_ids:
            selected_atom_nrs.append(atom.nr)
            selected_pdb_atoms.append(pdb_atom)

    summary = MoleculeTerminiSummary(
        itp_name="",
        total_atoms=len(atoms),
        candidate_atoms=sum(len(items) for items in candidates_by_chain.values()),
        chain_count=len(candidates_by_chain),
        restrained_atoms=len(selected_atom_nrs),
        min_coord=min_coord,
        max_coord=max_coord,
    )
    return summary, selected_atom_nrs, selected_pdb_atoms, selected_by_end


def write_summary(summary_path: Path, summaries: Sequence[MoleculeTerminiSummary]) -> None:
    ensure_parent(summary_path)
    with summary_path.open("w") as handle:
        handle.write(
            "itp\ttotal_atoms\tcandidate_atoms\tchain_count\trestrained_atoms\tmin_coord_A\tmax_coord_A\n"
        )
        for item in summaries:
            handle.write(
                f"{item.itp_name}\t{item.total_atoms}\t{item.candidate_atoms}\t"
                f"{item.chain_count}\t{item.restrained_atoms}\t"
                f"{item.min_coord:.3f}\t{item.max_coord:.3f}\n"
            )


def write_ndx_groups(index_path: Path, groups: Sequence[tuple[str, Sequence[int]]]) -> None:
    ensure_parent(index_path)
    with index_path.open("w") as handle:
        for group_name, atom_numbers in groups:
            handle.write(f"[ {group_name} ]\n")
            if atom_numbers:
                for start in range(0, len(atom_numbers), 15):
                    chunk = atom_numbers[start : start + 15]
                    handle.write(" ".join(str(atom_number) for atom_number in chunk) + "\n")
            handle.write("\n")


def parse_ndx_groups(index_path: Path) -> List[tuple[str, List[int]]]:
    groups: List[tuple[str, List[int]]] = []
    current_name: str | None = None
    current_atoms: List[int] = []

    with index_path.open() as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.startswith("[") and stripped.endswith("]"):
                if current_name is not None:
                    groups.append((current_name, current_atoms))
                current_name = stripped[1:-1].strip()
                current_atoms = []
                continue
            if current_name is None:
                continue
            current_atoms.extend(int(field) for field in stripped.split())

    if current_name is not None:
        groups.append((current_name, current_atoms))
    return groups


def merge_ndx_groups(
    index_path: Path,
    groups: Sequence[tuple[str, Sequence[int]]],
    purge_prefix: str | None = None,
) -> None:
    existing_groups = parse_ndx_groups(index_path) if index_path.exists() else []
    replacement_names = {name for name, _ in groups}

    merged_groups: List[tuple[str, Sequence[int]]] = []
    purge_stem = f"{purge_prefix}_COL" if purge_prefix else None

    for group_name, atom_numbers in existing_groups:
        if group_name in replacement_names:
            continue
        if purge_stem is not None and group_name.startswith(purge_stem):
            continue
        merged_groups.append((group_name, atom_numbers))

    for group_name, atom_numbers in groups:
        merged_groups.append((group_name, atom_numbers))

    write_ndx_groups(index_path, merged_groups)


def axis_to_pull_dim(axis: str) -> str:
    mapping = {
        "x": "Y N N",
        "y": "N Y N",
        "z": "N N Y",
    }
    return mapping[axis]


def write_pull_mdp_snippet(
    snippet_path: Path,
    pair_names: Sequence[tuple[str, str]],
    axis: str,
    pull_type: str,
    pull_k: str,
    pull_rate: str,
) -> None:
    ensure_parent(snippet_path)
    dim = axis_to_pull_dim(axis)
    lines: List[str] = [
        "; Auto-generated by add_termini_posres.py",
        "; Requires matching groups in the referenced .ndx file.",
        "pull                    = yes",
        f"pull-ngroups            = {len(pair_names) * 2}",
        f"pull-ncoords            = {len(pair_names)}",
        "pull-pbc-ref-prev-step-com = yes",
        "pull-nstxout            = 100",
        "pull-nstfout            = 100",
        "",
    ]

    for coord_index, (low_name, high_name) in enumerate(pair_names, start=1):
        low_group = (coord_index * 2) - 1
        high_group = coord_index * 2
        lines.extend(
            [
                f"pull-group{low_group}-name       = {low_name}",
                f"pull-group{high_group}-name       = {high_name}",
            ]
        )

    lines.append("")

    for coord_index, _pair in enumerate(pair_names, start=1):
        low_group = (coord_index * 2) - 1
        high_group = coord_index * 2
        lines.extend(
            [
                f"pull-coord{coord_index}-groups      = {low_group} {high_group}",
                f"pull-coord{coord_index}-geometry    = distance",
                f"pull-coord{coord_index}-dim         = {dim}",
                f"pull-coord{coord_index}-type        = {pull_type}",
            ]
        )
        if pull_type == "umbrella":
            lines.extend(
                [
                    f"pull-coord{coord_index}-start       = yes",
                    f"pull-coord{coord_index}-init        = 0.0",
                    f"pull-coord{coord_index}-rate        = {pull_rate}",
                    f"pull-coord{coord_index}-k           = {pull_k}",
                ]
            )
        else:
            lines.extend(
                [
                    f"pull-coord{coord_index}-k           = {pull_k}",
                ]
            )
        lines.append("")

    if pull_type == "constant-force":
        lines.extend(
            [
                "; Note: for extension along the chosen axis, pull-coord<n>-k is usually negative.",
                "",
            ]
        )

    snippet_path.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    if args.per_chain_count <= 0:
        raise ValueError("--per-chain-count must be positive")

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

    index_path = None
    if args.index_file:
        index_path = args.index_file
        if not index_path.is_absolute():
            index_path = topology_dir / index_path

    pull_mdp_snippet_path = None
    if args.pull_mdp_snippet:
        pull_mdp_snippet_path = args.pull_mdp_snippet
        if not pull_mdp_snippet_path.is_absolute():
            pull_mdp_snippet_path = topology_dir / pull_mdp_snippet_path

    if pull_mdp_snippet_path is not None and index_path is None:
        raise ValueError("--pull-mdp-snippet requires --index-file")
    if pull_mdp_snippet_path is not None and args.end != "both":
        raise ValueError("--pull-mdp-snippet currently requires --end both")
    if pull_mdp_snippet_path is not None and args.index_detail == "chain":
        raise ValueError(
            "--pull-mdp-snippet does not support --index-detail chain. "
            "Use aggregate, molecule, or all."
        )

    pdb_atoms = parse_pdb_atoms(pdb_path)
    parsed_itps = [(path, parse_itp_atoms(path)) for path in itp_paths]
    pdb_chunks = chunk_records(
        pdb_atoms,
        atom_counts=(len(atoms) for _, atoms in parsed_itps),
    )

    summaries: List[MoleculeTerminiSummary] = []
    all_restrained_pdb_atoms: List[PdbAtomRecord] = []
    all_low_end_atoms: List[int] = []
    all_high_end_atoms: List[int] = []
    ndx_groups: List[tuple[str, Sequence[int]]] = []
    molecule_pull_pairs: List[tuple[str, str]] = []

    for molecule_index, ((itp_path, atoms), pdb_chunk) in enumerate(zip(parsed_itps, pdb_chunks)):
        summary, restrained_atoms, restrained_pdb_atoms, selected_by_end = select_chain_termini(
            atoms=atoms,
            pdb_atoms=pdb_chunk,
            axis=args.axis,
            atom_name=args.atom_name,
            end_mode=args.end,
            per_chain_count=args.per_chain_count,
        )
        summary.itp_name = itp_path.name
        summaries.append(summary)
        all_restrained_pdb_atoms.extend(restrained_pdb_atoms)

        molecule_label = f"{args.index_prefix}_COL{molecule_index}"

        if args.end in {"both", "low"}:
            molecule_low: List[int] = []
            for chain_id in sorted(selected_by_end["low"]):
                chain_atoms = [atom.index + 1 for atom in selected_by_end["low"][chain_id]]
                molecule_low.extend(chain_atoms)
                if args.index_detail in {"chain", "all"}:
                    ndx_groups.append(
                        (f"{molecule_label}_{chain_id}_LOW", chain_atoms)
                    )
            all_low_end_atoms.extend(molecule_low)
            if args.index_detail in {"molecule", "all"}:
                ndx_groups.append((f"{molecule_label}_LOW", molecule_low))

        if args.end in {"both", "high"}:
            molecule_high: List[int] = []
            for chain_id in sorted(selected_by_end["high"]):
                chain_atoms = [atom.index + 1 for atom in selected_by_end["high"][chain_id]]
                molecule_high.extend(chain_atoms)
                if args.index_detail in {"chain", "all"}:
                    ndx_groups.append(
                        (f"{molecule_label}_{chain_id}_HIGH", chain_atoms)
                    )
            all_high_end_atoms.extend(molecule_high)
            if args.index_detail in {"molecule", "all"}:
                ndx_groups.append((f"{molecule_label}_HIGH", molecule_high))

        if args.end == "both" and args.index_detail in {"molecule", "all"}:
            molecule_pull_pairs.append((f"{molecule_label}_LOW", f"{molecule_label}_HIGH"))

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

    if index_path is not None:
        aggregate_groups: List[tuple[str, Sequence[int]]] = []
        if args.end in {"both", "low"}:
            aggregate_groups.append((f"{args.index_prefix}_LOW", all_low_end_atoms))
        if args.end in {"both", "high"}:
            aggregate_groups.append((f"{args.index_prefix}_HIGH", all_high_end_atoms))
        if args.end == "both":
            aggregate_groups.append(
                (
                    f"{args.index_prefix}_BOTH",
                    list(all_low_end_atoms) + list(all_high_end_atoms),
                )
            )
        exported_groups = aggregate_groups + ndx_groups
        if args.merge_index:
            merge_ndx_groups(
                index_path,
                exported_groups,
                purge_prefix=args.index_prefix,
            )
        else:
            write_ndx_groups(index_path, exported_groups)

    if pull_mdp_snippet_path is not None:
        if args.index_detail in {"molecule", "all"}:
            pair_names = molecule_pull_pairs
        else:
            pair_names = [(f"{args.index_prefix}_LOW", f"{args.index_prefix}_HIGH")]
        write_pull_mdp_snippet(
            snippet_path=pull_mdp_snippet_path,
            pair_names=pair_names,
            axis=args.axis,
            pull_type=args.pull_type,
            pull_k=args.pull_k,
            pull_rate=args.pull_rate,
        )

    action = "Would rewrite" if args.dry_run else "Rewrote"
    print(
        f"{action} {len(summaries)} ITP files in {topology_dir} using axis={args.axis}, "
        f"per-chain-count={args.per_chain_count}, end={args.end}, atom={args.atom_name}."
    )
    for item in summaries:
        print(
            f"{item.itp_name}: restrained {item.restrained_atoms}/{item.candidate_atoms} "
            f"candidate atoms across {item.chain_count} chains between "
            f"{item.min_coord:.3f} A and {item.max_coord:.3f} A"
        )
    if index_path is not None:
        action_label = "Merged" if args.merge_index else "Wrote"
        print(f"{action_label} GROMACS index groups to {index_path}")
    if pull_mdp_snippet_path is not None:
        print(f"Wrote pull-code .mdp snippet to {pull_mdp_snippet_path}")


if __name__ == "__main__":
    main()
