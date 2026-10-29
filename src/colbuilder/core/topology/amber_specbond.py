"""Opt-in Amber topology generation using GROMACS special-bond templates.

This backend deliberately keeps :mod:`colbuilder.core.topology.amber` intact.
It splits the final fibril in its on-disk order, lets ``pdb2gmx`` generate the
entire bonded topology for each connected component, and verifies that every
crosslink expected from the input geometry is present in that topology.
"""

from __future__ import annotations

import asyncio
import math
import os
import re
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from colbuilder.core.geometry.crosslink import Crosslink
from colbuilder.core.geometry.system import System
from colbuilder.core.topology.amber import Amber
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.exceptions import TopologyGenerationError
from colbuilder.core.utils.files import FileManager
from colbuilder.core.utils.logger import setup_logger


LOG = setup_logger(__name__)

FF = "amber99sb-star-ildnp"
FF_DIR = f"{FF}.ff"
GROUPING_CUTOFF_A = 5.0
TOPOLOGY_CUTOFF_A = 5.0

SECTION_RE = re.compile(r"^\s*\[\s*([^\]]+?)\s*\]")
INCLUDE_RE = re.compile(r'^(\s*#include\s+")([^"]+)(".*)$')

XLINK_ATOMS: Dict[str, Set[str]] = {
    "LYX": {"C13", "C12"},
    "LXY": {"C13", "C12"},
    "LYY": {"C13", "C12"},
    "LXX": {"C13", "C12"},
    "LY3": {"CG"},
    "LX3": {"CG"},
    "L3Y": {"CG"},
    "L2Y": {"CG"},
    "L3X": {"CG"},
    "L2X": {"CG"},
    "LY2": {"CB"},
    "LX2": {"CB"},
    "L4Y": {"CE"},
    "L4X": {"CE"},
    "LY4": {"CE"},
    "LX4": {"CE"},
    "L5Y": {"NZ"},
    "L5X": {"NZ"},
    "LY5": {"NZ"},
    "LX5": {"NZ"},
    "LGX": {"CE"},
    "LPS": {"CE"},
    "AGS": {"NZ"},
    "APD": {"NZ"},
    "LZD": {"CE"},
    "LZS": {"NZ1"},
}

TRIVALENT_RESIDUES = {
    "LYX", "LXY", "LYY", "LXX", "LY3", "LX3", "L3Y", "L2Y",
    "L3X", "L2X", "LY2", "LX2",
}
CENTRAL_RESIDUES = {"LYX", "LXY", "LYY", "LXX"}

Marker = Dict[str, object]
Hit = Dict[str, object]


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, first: int, second: int) -> None:
        first_root, second_root = self.find(first), self.find(second)
        if first_root != second_root:
            self.parent[second_root] = first_root

    def groups(self) -> List[List[int]]:
        grouped: Dict[int, List[int]] = {}
        for item in range(len(self.parent)):
            grouped.setdefault(self.find(item), []).append(item)
        return sorted((sorted(group) for group in grouped.values()), key=lambda g: g[0])


def _record_type(line: str) -> str:
    return line[:6].strip().upper()


def _section(line: str) -> Optional[str]:
    match = SECTION_RE.match(line)
    return match.group(1).strip().lower() if match else None


def split_final_fibril(
    pdb: Path,
    chains_per_collagen: int = 3,
) -> Tuple[List[str], List[List[List[str]]]]:
    """Split a final fibril by chain-ID changes and assign sequential IDs."""
    header: List[str] = []
    chains: List[List[str]] = []
    current: List[str] = []
    current_chain: Optional[str] = None
    seen_atoms = False

    def flush() -> None:
        nonlocal current, current_chain
        if current:
            current = [line for line in current if _record_type(line) != "TER"]
            current.append("TER\n")
            chains.append(current)
        current = []
        current_chain = None

    with pdb.open() as handle:
        for line in handle:
            record = _record_type(line)
            if record in {"ATOM", "HETATM"}:
                seen_atoms = True
                chain_id = line[21:22]
                if current and current_chain is not None and chain_id != current_chain:
                    flush()
                current.append(line)
                current_chain = chain_id
            elif record == "ANISOU" and current:
                current.append(line)
            elif not seen_atoms and record in {
                "HEADER", "TITLE", "REMARK", "CRYST1", "SCALE1", "SCALE2", "SCALE3"
            }:
                header.append(line)
    flush()

    if not chains:
        raise RuntimeError(f"No coordinate chains found in final fibril {pdb}")
    if len(chains) % chains_per_collagen:
        raise RuntimeError(
            f"Final fibril has {len(chains)} chemical chains, which is not divisible "
            f"by {chains_per_collagen}; check chain IDs in {pdb}"
        )
    collagens = [
        chains[index:index + chains_per_collagen]
        for index in range(0, len(chains), chains_per_collagen)
    ]
    return header, collagens


def _markers(collagens: Sequence[Sequence[Sequence[str]]]) -> List[Marker]:
    markers: List[Marker] = []
    for collagen_id, collagen in enumerate(collagens):
        for chain_block in collagen:
            for line in chain_block:
                if _record_type(line) not in {"ATOM", "HETATM"}:
                    continue
                resname = line[17:20].strip()
                atom = line[12:16].strip()
                if atom not in XLINK_ATOMS.get(resname, set()):
                    continue
                markers.append(
                    {
                        "collagen": collagen_id,
                        "resname": resname,
                        "atom": atom,
                        "resid": line[22:26].strip(),
                        "chain": line[21:22].strip(),
                        "xyz": (
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54]),
                        ),
                    }
                )
    return markers


def _as_crosslink(marker: Marker) -> Crosslink:
    resname = str(marker["resname"])
    return Crosslink(
        resid=str(marker["resid"]),
        resname=resname,
        chain=str(marker["chain"]),
        position=list(marker["xyz"]),
        type="T" if resname in TRIVALENT_RESIDUES else "D",
        model_id=int(marker["collagen"]),
        atom=str(marker["atom"]),
    )


def _chemistry(marker_a: Marker, marker_b: Marker) -> Tuple[str, str, str, str]:
    return (
        str(marker_a["resname"]),
        str(marker_a["atom"]),
        str(marker_b["resname"]),
        str(marker_b["atom"]),
    )


def _distance(marker_a: Marker, marker_b: Marker) -> float:
    xyz_a = marker_a["xyz"]
    xyz_b = marker_b["xyz"]
    return math.sqrt(sum((xyz_a[i] - xyz_b[i]) ** 2 for i in range(3)))


def _compatible_hits(markers: Sequence[Marker], amber: Amber) -> List[Hit]:
    hits: List[Hit] = []
    for index, first in enumerate(markers):
        for second in markers[index + 1:]:
            if first["collagen"] == second["collagen"]:
                continue
            if not amber._are_compatible_crosslinks(  # reuse legacy chemistry rules
                _as_crosslink(first), _as_crosslink(second)
            ):
                continue
            distance = _distance(first, second)
            if distance < TOPOLOGY_CUTOFF_A:
                hits.append({"a": first, "b": second, "distance_A": distance})
    return hits


def _marker_identity(marker: Marker) -> Tuple[object, ...]:
    return (
        marker["collagen"], marker["resname"], marker["resid"],
        marker["chain"], marker["atom"],
    )


def _hit_description(hit: Hit) -> str:
    first, second = hit["a"], hit["b"]
    return (
        f"col_{first['collagen']}:{first['resname']}{first['resid']}"
        f"{first['chain']}/{first['atom']} <-> "
        f"col_{second['collagen']}:{second['resname']}{second['resid']}"
        f"{second['chain']}/{second['atom']} {hit['distance_A']:.3f} A"
    )


def _validate_unambiguous(hits: Sequence[Hit]) -> None:
    partners: Dict[Tuple[object, ...], Set[Tuple[object, ...]]] = {}
    for hit in hits:
        first = _marker_identity(hit["a"])
        second = _marker_identity(hit["b"])
        partners.setdefault(first, set()).add(second)
        partners.setdefault(second, set()).add(first)
    ambiguous = [identity for identity, values in partners.items() if len(values) > 1]
    if ambiguous:
        raise RuntimeError(
            "Ambiguous crosslink chemistry: reactive atom(s) have multiple partners "
            f"within {TOPOLOGY_CUTOFF_A:g} A: {ambiguous}"
        )


def _pyd_arm(hit: Hit) -> Optional[str]:
    for central, arm in ((hit["a"], hit["b"]), (hit["b"], hit["a"])):
        if central["resname"] not in CENTRAL_RESIDUES:
            continue
        # The arm residue for each ring carbon is checked by
        # Amber._are_compatible_crosslinks; PYL/DPL bond C12 via CG, not CB.
        if central["atom"] in ("C13", "C12"):
            return str(central["atom"])
    return None


def _validate_complete_pyd(hits: Sequence[Hit]) -> int:
    centres: Dict[Tuple[object, ...], Set[str]] = {}
    for hit in hits:
        arm = _pyd_arm(hit)
        if arm is None:
            continue
        central = next(
            marker for marker in (hit["a"], hit["b"])
            if marker["resname"] in CENTRAL_RESIDUES
        )
        identity = (
            central["collagen"], central["resname"], central["resid"], central["chain"]
        )
        centres.setdefault(identity, set()).add(arm)
    incomplete = {
        identity: sorted({"C13", "C12"} - arms)
        for identity, arms in centres.items()
        if arms != {"C13", "C12"}
    }
    if incomplete:
        raise RuntimeError(f"Incomplete trivalent centre(s); missing ring bond: {incomplete}")
    return len(centres)


def read_specbond_records(path: Path) -> Tuple[List[List[str]], List[str]]:
    records: List[List[str]] = []
    comments: List[str] = []
    with path.open() as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped:
                continue
            if stripped.startswith((";", "#")):
                comments.append(raw.rstrip("\n"))
                continue
            fields = raw.split(";", 1)[0].split()
            if len(fields) == 1 and fields[0].isdigit():
                continue
            if len(fields) >= 9:
                records.append(fields[:9])
    return records, comments


def _record_matches(
    record: Sequence[str], chemistry: Tuple[str, str, str, str]
) -> bool:
    forward = (record[0], record[1], record[3], record[4])
    reverse = (record[3], record[4], record[0], record[1])
    return chemistry in {forward, reverse}


def recognition_reference_clusters_nm(
    distances_nm: Iterable[float],
) -> List[Tuple[float, List[float]]]:
    """Cluster distances into the fewest practical GROMACS +/-10% windows."""
    values = sorted(float(value) for value in distances_nm)
    if not values:
        return []
    clusters: List[Tuple[float, List[float]]] = []
    current: List[float] = []
    low = high = 0.0
    for value in values:
        value_low, value_high = value / 1.1, value / 0.9
        if not current:
            current = [value]
            low, high = value_low, value_high
            continue
        new_low, new_high = max(low, value_low), min(high, value_high)
        if new_low <= new_high:
            current.append(value)
            low, high = new_low, new_high
        else:
            clusters.append(((low + high) / 2.0, current))
            current = [value]
            low, high = value_low, value_high
    clusters.append(((low + high) / 2.0, current))
    return clusters


def write_component_specbond(
    base_specbond: Path,
    output: Path,
    component_hits: Sequence[Hit],
) -> None:
    """Write component-specific recognition windows from observed distances."""
    records, comments = read_specbond_records(base_specbond)
    by_chemistry: Dict[Tuple[Tuple[str, str], Tuple[str, str]], List[float]] = {}
    for hit in component_hits:
        first, second = hit["a"], hit["b"]
        key = tuple(sorted((
            (str(first["resname"]), str(first["atom"])),
            (str(second["resname"]), str(second["atom"])),
        )))
        by_chemistry.setdefault(key, []).append(float(hit["distance_A"]) / 10.0)

    for key, distances in by_chemistry.items():
        chemistry = (key[0][0], key[0][1], key[1][0], key[1][1])
        matching = [
            index for index, record in enumerate(records)
            if _record_matches(record, chemistry)
        ]
        if not matching:
            raise RuntimeError(
                "No specbond.dat chemistry for required crosslink "
                f"{chemistry[0]}/{chemistry[1]} <-> {chemistry[2]}/{chemistry[3]}"
            )
        insertion = matching[0]
        template = list(records[insertion])
        for index in reversed(matching):
            del records[index]
        for offset, (reference, _) in enumerate(
            recognition_reference_clusters_nm(distances)
        ):
            replacement = list(template)
            replacement[6] = f"{reference:.4f}"
            records.insert(insertion + offset, replacement)

    with output.open("w") as handle:
        handle.write(f"{len(records)}\n")
        for record in records:
            handle.write("  ".join(record) + "\n")
        if comments:
            handle.write("\n" + "\n".join(comments) + "\n")


def _chemistries(records: Sequence[Sequence[str]]) -> Set[Tuple[str, str, str, str]]:
    allowed: Set[Tuple[str, str, str, str]] = set()
    for record in records:
        forward = (record[0], record[1], record[3], record[4])
        allowed.add(forward)
        allowed.add((forward[2], forward[3], forward[0], forward[1]))
    return allowed


def build_sequential_components(
    collagens: Sequence[Sequence[Sequence[str]]],
    amber: Amber,
    specbond: Path,
) -> Tuple[List[List[int]], List[Hit], int]:
    """Build sequential component IDs and validate their expected chemistry."""
    records, _ = read_specbond_records(specbond)
    allowed = _chemistries(records)
    hits = _compatible_hits(_markers(collagens), amber)

    missing_chemistry = [
        hit for hit in hits
        if _chemistry(hit["a"], hit["b"]) not in allowed
    ]
    if missing_chemistry:
        raise RuntimeError(
            "Required crosslink chemistry is absent from specbond.dat: "
            + "; ".join(_hit_description(hit) for hit in missing_chemistry)
        )

    _validate_unambiguous(hits)
    complete_pyd = _validate_complete_pyd(hits)

    union = UnionFind(len(collagens))
    for hit in hits:
        if float(hit["distance_A"]) < GROUPING_CUTOFF_A:
            union.union(int(hit["a"]["collagen"]), int(hit["b"]["collagen"]))
    groups = union.groups()
    component = {
        collagen_id: group_index
        for group_index, group in enumerate(groups)
        for collagen_id in group
    }
    missing_components = [
        hit for hit in hits
        if component[int(hit["a"]["collagen"])]
        != component[int(hit["b"]["collagen"])]
    ]
    if missing_components:
        raise RuntimeError(
            "Crosslink component is missing: compatible atom pair(s) within the "
            f"{TOPOLOGY_CUTOFF_A:g} A topology cutoff were not connected by the "
            f"{GROUPING_CUTOFF_A:g} A component cutoff: "
            + "; ".join(_hit_description(hit) for hit in missing_components)
        )
    return groups, hits, complete_pyd


def _write_component_pdb(
    path: Path,
    header: Sequence[str],
    collagens: Sequence[Sequence[Sequence[str]]],
    ids: Sequence[int],
) -> None:
    with path.open("w") as handle:
        handle.writelines(header)
        for collagen_id in ids:
            for chain in collagens[collagen_id]:
                handle.writelines(chain)
        handle.write("END\n")


def _rename_moleculetype(text: str, name: str) -> str:
    lines = text.splitlines(keepends=True)
    for index, line in enumerate(lines):
        if _section(line) != "moleculetype":
            continue
        for data_index in range(index + 1, len(lines)):
            if _section(lines[data_index]):
                break
            stripped = lines[data_index].strip()
            if not stripped or stripped.startswith((";", "#")):
                continue
            fields = stripped.split()
            newline = "\n" if lines[data_index].endswith("\n") else ""
            lines[data_index] = f"{name}  {fields[1]}{newline}"
            return "".join(lines)
    raise RuntimeError("pdb2gmx output has no [ moleculetype ] section")


def _inline_molecule(topology: Path) -> Optional[str]:
    lines = topology.read_text().splitlines(keepends=True)
    start: Optional[int] = None
    for index, line in enumerate(lines):
        if _section(line) == "moleculetype":
            start = index
            break
    if start is None:
        return None
    end = len(lines)
    for index in range(start + 1, len(lines)):
        lower = lines[index].lower()
        if (
            _section(lines[index]) in {"system", "molecules"}
            or "include water topology" in lower
            or "include topology for ions" in lower
        ):
            end = index
            break
    return "".join(lines[start:end]).rstrip() + "\n"


def _moleculetype_name(text: str) -> Optional[str]:
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if _section(line) != "moleculetype":
            continue
        for data in lines[index + 1:]:
            if _section(data):
                break
            stripped = data.strip()
            if stripped and not stripped.startswith((";", "#")):
                return stripped.split()[0]
    return None


def _normalize_posre_include(text: str, posre_name: str) -> str:
    normalized: List[str] = []
    for line in text.splitlines(keepends=True):
        match = INCLUDE_RE.match(line)
        if match and "posre" in Path(match.group(2)).name.lower():
            line = match.group(1) + posre_name + match.group(3)
        normalized.append(line)
    return "".join(normalized)


def keep_generated_molecule_topology(
    generated_top: Path,
    workdir: Path,
    output_itp: Path,
    molecule_name: str,
    posre_name: str,
) -> None:
    """Keep pdb2gmx's molecule topology without synthesizing bonded terms."""
    candidates: List[Path] = []
    for candidate in workdir.glob("*.itp"):
        try:
            if _moleculetype_name(candidate.read_text()):
                candidates.append(candidate)
        except UnicodeDecodeError:
            continue
    if len(candidates) == 1:
        text = candidates[0].read_text()
    elif not candidates:
        text = _inline_molecule(generated_top)
        if text is None:
            raise RuntimeError("pdb2gmx generated no molecule topology")
    else:
        raise RuntimeError(
            f"pdb2gmx generated multiple molecule ITPs in {workdir}: {candidates}"
        )
    text = _rename_moleculetype(text, molecule_name)
    output_itp.write_text(_normalize_posre_include(text, posre_name))


def _parse_gro_atoms(gro: Path) -> List[Marker]:
    lines = gro.read_text().splitlines()
    count = int(lines[1].strip())
    atoms: List[Marker] = []
    for line in lines[2:2 + count]:
        atoms.append(
            {
                "resname": line[5:10].strip(),
                "atom": line[10:15].strip(),
                "index": int(line[15:20].strip()),
                "xyz": (
                    float(line[20:28]) * 10.0,
                    float(line[28:36]) * 10.0,
                    float(line[36:44]) * 10.0,
                ),
            }
        )
    return atoms


def _nearest_atom(atoms: Sequence[Marker], target: Marker) -> int:
    best_index: Optional[int] = None
    best_distance = float("inf")
    for atom in atoms:
        if atom["resname"] != target["resname"] or atom["atom"] != target["atom"]:
            continue
        distance = _distance(atom, target)
        if distance < best_distance:
            best_index, best_distance = int(atom["index"]), distance
    if best_index is None or best_distance > 0.25:
        raise RuntimeError(
            "Could not map expected crosslink atom "
            f"{target['resname']}{target['resid']}{target['chain']}/{target['atom']} "
            f"to pdb2gmx coordinates within 0.25 A (nearest {best_distance:.3f} A)"
        )
    return best_index


def _parse_itp_bonds(itp: Path) -> Set[Tuple[int, int]]:
    bonds: Set[Tuple[int, int]] = set()
    in_bonds = False
    with itp.open() as handle:
        for raw in handle:
            section = _section(raw)
            if section is not None:
                in_bonds = section == "bonds"
                continue
            if not in_bonds:
                continue
            fields = raw.split(";", 1)[0].split()
            if len(fields) < 2:
                continue
            try:
                first, second = int(fields[0]), int(fields[1])
            except ValueError:
                continue
            bonds.add(tuple(sorted((first, second))))
    return bonds


def verify_component_crosslinks(
    gro: Path,
    itp: Path,
    component_hits: Sequence[Hit],
) -> int:
    """Fail unless every expected crosslink is in pdb2gmx's [ bonds ]."""
    if not component_hits:
        return 0
    atoms = _parse_gro_atoms(gro)
    bonds = _parse_itp_bonds(itp)
    for hit in component_hits:
        first = _nearest_atom(atoms, hit["a"])
        second = _nearest_atom(atoms, hit["b"])
        if tuple(sorted((first, second))) not in bonds:
            raise RuntimeError(
                "pdb2gmx did not create expected crosslink bond: "
                f"{_hit_description(hit)} -> ITP atoms {first}-{second}"
            )
    return len(component_hits)


async def _run_pdb2gmx(
    pdb: Path,
    top: Path,
    gro: Path,
    posre: Path,
    workdir: Path,
    gmxlib: Path,
    water_type: str,
) -> None:
    command = [
        "gmx", "pdb2gmx", "-f", str(pdb), "-o", str(gro), "-p", str(top),
        "-i", str(posre), "-ff", FF, "-water", water_type, "-chainsep", "id",
        "-merge", "all", "-ignh",
    ]
    environment = os.environ.copy()
    environment["GMXLIB"] = str(gmxlib)
    process = await asyncio.create_subprocess_exec(
        *command,
        cwd=str(workdir),
        env=environment,
        stdin=asyncio.subprocess.DEVNULL,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
    )
    output, _ = await process.communicate()
    log = workdir / "pdb2gmx.log"
    log.write_bytes(output or b"")
    if process.returncode:
        raise RuntimeError(
            f"GROMACS pdb2gmx failed for component {workdir.name}; see {log}"
        )


def _component_gmxlib(
    workdir: Path,
    base_gmxlib: Path,
    base_specbond: Path,
    component_hits: Sequence[Hit],
) -> Path:
    if not component_hits:
        return base_gmxlib
    runtime = workdir / "gmxlib"
    runtime.mkdir()
    os.symlink(base_gmxlib / FF_DIR, runtime / FF_DIR, target_is_directory=True)
    os.symlink(base_gmxlib / "residuetypes.dat", runtime / "residuetypes.dat")
    write_component_specbond(base_specbond, runtime / "specbond.dat", component_hits)
    return runtime


@timeit
async def build_amber99_specbond(
    system: System,
    config: ColbuilderConfig,
    file_manager: Optional[FileManager] = None,
) -> Amber:
    """Build an Amber topology entirely from pdb2gmx special-bond output."""
    del file_manager  # Kept for the same backend signature as build_amber99.
    working_dir = Path.cwd()
    fibril = Path(system.pdb_fibril).resolve()
    source_ff = config.FORCE_FIELD_DIR / FF_DIR
    source_specbond = source_ff / "specbond.dat"
    source_residuetypes = source_ff / "residuetypes.dat"
    amber = Amber(system=system, ff=FF, water_type=config.water_type)
    workroot = working_dir / ".specbond_work"

    try:
        if not fibril.is_file():
            raise RuntimeError(
                "Specbond topology mode requires the final fibril PDB; "
                f"not found: {fibril}"
            )
        for resource in (source_ff, source_specbond, source_residuetypes):
            if not resource.exists():
                raise RuntimeError(f"Required Amber specbond resource is missing: {resource}")

        copied_ff = working_dir / FF_DIR
        if not copied_ff.exists():
            shutil.copytree(source_ff, copied_ff)
        shutil.copy2(source_specbond, working_dir / "specbond.dat")
        shutil.copy2(source_residuetypes, working_dir / "residuetypes.dat")
        base_specbond = working_dir / "specbond.dat"

        header, collagens = split_final_fibril(fibril)
        groups, hits, complete_pyd = build_sequential_components(
            collagens, amber, base_specbond
        )
        connected = sum(1 for group in groups if len(group) > 1)
        LOG.info(
            f"Specbond topology: {len(collagens)} sequential collagens, "
            f"{len(groups)} components ({connected} crosslink-connected), "
            f"{len(hits)} expected bonds"
        )

        if workroot.exists():
            shutil.rmtree(workroot)
        workroot.mkdir()
        processed_groups: List[Tuple[str, str]] = []
        verified_total = 0

        for group_index, ids in enumerate(groups):
            group_id = "_".join(str(item) for item in ids)
            molecule_name = f"col_{group_id}"
            component_dir = workroot / f"{group_index:03d}_{molecule_name}"
            component_dir.mkdir()
            component_pdb = component_dir / "component.pdb"
            generated_top = component_dir / "component.top"
            generated_gro = component_dir / "component.gro"
            posre = working_dir / f"posre_{group_id}.itp"
            final_itp = working_dir / f"col_{group_id}.itp"
            final_gro = working_dir / f"col_{group_id}.gro"
            component_hits = [
                hit for hit in hits
                if int(hit["a"]["collagen"]) in ids
                and int(hit["b"]["collagen"]) in ids
            ]

            _write_component_pdb(component_pdb, header, collagens, ids)
            runtime = _component_gmxlib(
                component_dir, working_dir, base_specbond, component_hits
            )
            await _run_pdb2gmx(
                component_pdb.resolve(), generated_top.resolve(),
                generated_gro.resolve(), posre.resolve(),
                component_dir.resolve(), runtime.resolve(), config.water_type,
            )
            keep_generated_molecule_topology(
                generated_top, component_dir, final_itp, molecule_name, posre.name
            )
            shutil.copy2(generated_gro, final_gro)
            amber.ensure_posre_include(final_itp, group_id)
            verified_total += verify_component_crosslinks(
                final_gro, final_itp, component_hits
            )
            processed_groups.append(("specbond", group_id))

        if verified_total != len(hits):
            raise RuntimeError(
                f"Crosslink verification incomplete: {verified_total}/{len(hits)} bonds"
            )

        amber.write_topology(
            str(working_dir / f"collagen_fibril_{config.species}.top"),
            processed_groups,
        )
        amber.write_gro(
            str(working_dir / f"collagen_fibril_{config.species}.gro"),
            processed_groups,
        )

        amber.component_count = len(groups)
        amber.connected_component_count = connected
        amber.verified_bond_count = verified_total
        amber.complete_pyd_centres = complete_pyd
        LOG.info(
            f"Specbond topology complete: {verified_total} verified crosslink bonds; "
            f"{complete_pyd} complete trivalent centres"
        )
        if not config.topology_debug:
            shutil.rmtree(workroot)
        return amber

    except TopologyGenerationError:
        raise
    except Exception as error:
        raise TopologyGenerationError(
            message=f"Amber specbond topology generation failed: {error}",
            original_error=error,
            error_code="TOP_ERR_001",
            context={"fibril": str(fibril), "mode": "specbond"},
        ) from error
