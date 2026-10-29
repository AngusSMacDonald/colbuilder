import csv
from pathlib import Path
import sys
import types

import pytest

# MODELLER is an optional proprietary sequence dependency and is not needed by
# these topology/CLI parsing tests. Stub its import surface before loading the
# top-level CLI module.
if "modeller" not in sys.modules:
    modeller = types.ModuleType("modeller")
    modeller.Environ = object
    automodel = types.ModuleType("modeller.automodel")
    automodel.AutoModel = object
    scripts = types.ModuleType("modeller.scripts")
    scripts.complete_pdb = lambda *args, **kwargs: None
    sys.modules.update(
        {
            "modeller": modeller,
            "modeller.automodel": automodel,
            "modeller.scripts": scripts,
        }
    )

from colbuilder.colbuilder import main
from colbuilder.core.geometry.crosslink import Crosslink
from colbuilder.core.topology.amber import Amber
from colbuilder.core.topology.amber_specbond import (
    XLINK_ATOMS,
    _validate_complete_pyd,
    count_crosslink_bonds,
    read_specbond_records,
    verify_component_crosslinks,
    write_component_specbond,
)
from colbuilder.core.utils.config import ColbuilderConfig, load_yaml_config


def _marker(collagen, resname, atom, resid, chain, xyz):
    return {
        "collagen": collagen,
        "resname": resname,
        "atom": atom,
        "resid": str(resid),
        "chain": chain,
        "xyz": xyz,
    }


def _hit(first, second, distance):
    return {"a": first, "b": second, "distance_A": distance}


def _base_specbond(path: Path) -> Path:
    path.write_text(
        "3\n"
        "L4Y CE 1 L5Y NZ 1 0.152 L4Y L5Y\n"
        "LY2 CB 1 LYX C12 1 0.152 LY2 LYX\n"
        "LY3 CG 1 LYX C13 1 0.152 LY3 LYX\n"
    )
    return path


def _reference_lengths(path: Path, chemistry):
    records, _ = read_specbond_records(path)
    output = []
    for record in records:
        forward = (record[0], record[1], record[3], record[4])
        reverse = (record[3], record[4], record[0], record[1])
        if chemistry in (forward, reverse):
            output.append(float(record[6]))
    return output


def test_amber_topology_mode_defaults_to_legacy():
    config = ColbuilderConfig(species="homo_sapiens")
    assert config.amber_topology_mode == "legacy"
    assert config.water_type == "tip3p"


def test_yaml_and_cli_select_specbond(tmp_path):
    yaml_path = tmp_path / "config.yaml"
    yaml_path.write_text(
        "species: homo_sapiens\n"
        "force_field: amber99\n"
        "amber_topology_mode: specbond\n"
        "water_type: tip4p\n"
    )
    config = ColbuilderConfig(**load_yaml_config(yaml_path))
    assert config.amber_topology_mode == "specbond"
    assert config.water_type == "tip4p"

    context = main.make_context(
        "colbuilder",
        [
            "--species", "homo_sapiens",
            "--force_field", "amber99",
            "--amber-topology-mode", "specbond",
        ],
    )
    assert context.params["amber_topology_mode"] == "specbond"


def test_specbond_mode_rejects_martini():
    with pytest.raises(Exception, match="only supported.*amber99"):
        ColbuilderConfig(
            species="homo_sapiens",
            topology_generator=True,
            force_field="martini3",
            amber_topology_mode="specbond",
        )


def test_hlknl_recognition_windows_follow_observed_distances(tmp_path):
    base = _base_specbond(tmp_path / "base.dat")
    output = tmp_path / "component.dat"
    aldehyde = _marker(0, "L4Y", "CE", 947, "A", (0.0, 0.0, 0.0))
    amine_a = _marker(1, "L5Y", "NZ", 9, "C", (1.81, 0.0, 0.0))
    amine_b = _marker(2, "L5Y", "NZ", 1047, "C", (2.41, 0.0, 0.0))
    write_component_specbond(
        base,
        output,
        [_hit(aldehyde, amine_a, 1.81), _hit(aldehyde, amine_b, 2.41)],
    )

    references = _reference_lengths(output, ("L4Y", "CE", "L5Y", "NZ"))
    assert len(references) == 2
    assert any(0.9 * ref <= 0.181 <= 1.1 * ref for ref in references)
    assert any(0.9 * ref <= 0.241 <= 1.1 * ref for ref in references)


def test_pyd_recognition_windows_cover_both_arms(tmp_path):
    base = _base_specbond(tmp_path / "base.dat")
    output = tmp_path / "component.dat"
    c12 = _marker(0, "LYX", "C12", 947, "A", (0.0, 0.0, 0.0))
    c13 = _marker(0, "LYX", "C13", 947, "A", (0.0, 1.0, 0.0))
    cb = _marker(1, "LY2", "CB", 9, "C", (2.35, 0.0, 0.0))
    cg = _marker(1, "LY3", "CG", 10, "C", (1.75, 1.0, 0.0))
    write_component_specbond(
        base,
        output,
        [_hit(c12, cb, 2.35), _hit(c13, cg, 1.75)],
    )

    c12_refs = _reference_lengths(output, ("LYX", "C12", "LY2", "CB"))
    c13_refs = _reference_lengths(output, ("LYX", "C13", "LY3", "CG"))
    assert len(c12_refs) == len(c13_refs) == 1
    assert 0.9 * c12_refs[0] <= 0.235 <= 1.1 * c12_refs[0]
    assert 0.9 * c13_refs[0] <= 0.175 <= 1.1 * c13_refs[0]


@pytest.mark.parametrize(
    "central, arm12, atom12, arm13",
    [("LYX", "LY2", "CB", "LY3"), ("LXX", "LX2", "CB", "LX3"),
     ("LXY", "L2Y", "CG", "L3Y"), ("LYY", "L2X", "CG", "L3X")],
)
def test_trivalent_centres_are_bundled_and_validated(central, arm12, atom12, arm13):
    bundled = (
        Path(__file__).parents[1]
        / "src/colbuilder/data/topology/amber99sb-star-ildnp.ff/specbond.dat"
    )
    records, _ = read_specbond_records(bundled)
    chemistries = {(r[0], r[1], r[3], r[4]) for r in records}
    assert (arm12, atom12, central, "C12") in chemistries
    assert (arm13, "CG", central, "C13") in chemistries

    c12 = _marker(0, central, "C12", 947, "A", (0.0, 0.0, 0.0))
    c13 = _marker(0, central, "C13", 947, "A", (0.0, 1.0, 0.0))
    first = _marker(1, arm12, atom12, 9, "C", (2.5, 0.0, 0.0))
    second = _marker(2, arm13, "CG", 5, "B", (2.5, 1.0, 0.0))
    assert _validate_complete_pyd([_hit(c12, first, 2.5), _hit(c13, second, 2.5)]) == 1
    with pytest.raises(RuntimeError, match="Incomplete trivalent"):
        _validate_complete_pyd([_hit(c12, first, 2.5)])


def _gro_atom(resid, resname, atom, index, x, y, z):
    return (
        f"{resid:5d}{resname:<5}{atom:>5}{index:5d}"
        f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
    )


def test_bond_verification_succeeds_and_missing_bond_fails(tmp_path):
    gro = tmp_path / "component.gro"
    gro.write_text(
        "test\n2\n"
        + _gro_atom(1, "L4Y", "CE", 1, 0.000, 0.000, 0.000)
        + _gro_atom(2, "L5Y", "NZ", 2, 0.181, 0.000, 0.000)
        + "1.0 1.0 1.0\n"
    )
    itp = tmp_path / "component.itp"
    atoms = (
        "[ atoms ]\n"
        "1 CT 1 L4Y CE 1 0.0 12.01\n"
        "2 N3 2 L5Y NZ 2 0.0 14.01\n"
    )
    itp.write_text(atoms + "[ bonds ]\n1 2 1\n")
    hit = _hit(
        _marker(0, "L4Y", "CE", 947, "A", (0.0, 0.0, 0.0)),
        _marker(1, "L5Y", "NZ", 9, "C", (1.81, 0.0, 0.0)),
        1.81,
    )
    assert verify_component_crosslinks(gro, itp, [hit]) == 1

    itp.write_text(atoms + "[ bonds ]\n")
    with pytest.raises(RuntimeError, match="did not create expected crosslink bond"):
        verify_component_crosslinks(gro, itp, [hit])


def test_gro_atom_index_uses_line_position_past_99999(tmp_path):
    count = 100002
    gro = tmp_path / "component.gro"
    gro.write_text(
        "test\n"
        f"{count}\n"
        + "".join(
            _gro_atom(i % 100000, "GLY", "CA", i % 100000, 9.0, 9.0, 9.0)
            for i in range(1, count - 1)
        )
        + _gro_atom(100001 % 100000, "L4Y", "CE", 100001 % 100000, 0.0, 0.0, 0.0)
        + _gro_atom(100002 % 100000, "L5Y", "NZ", 100002 % 100000, 0.181, 0.0, 0.0)
        + "1.0 1.0 1.0\n"
    )
    itp = tmp_path / "component.itp"
    atoms = "".join(
        f"{i} CT {i} GLY CA {i} 0.0 12.01\n" for i in range(1, count - 1)
    )
    atoms += f"{count - 1} CT {count - 1} L4Y CE {count - 1} 0.0 12.01\n"
    atoms += f"{count} N3 {count} L5Y NZ {count} 0.0 14.01\n"
    itp.write_text("[ atoms ]\n" + atoms + "[ bonds ]\n100001 100002 1\n")
    hit = _hit(
        _marker(0, "L4Y", "CE", 947, "A", (0.0, 0.0, 0.0)),
        _marker(1, "L5Y", "NZ", 9, "C", (1.81, 0.0, 0.0)),
        1.81,
    )
    assert verify_component_crosslinks(gro, itp, [hit]) == 1


def test_specbond_reference_covers_window_edges(tmp_path):
    base = _base_specbond(tmp_path / "base.dat")
    output = tmp_path / "component.dat"
    aldehyde = _marker(0, "L4Y", "CE", 947, "A", (0.0, 0.0, 0.0))
    amine_a = _marker(1, "L5Y", "NZ", 9, "C", (1.80, 0.0, 0.0))
    amine_b = _marker(2, "L5Y", "NZ", 1047, "C", (1.80 * 1.1 / 0.9 * 0.99999, 0.0, 0.0))
    distances = [1.80, 1.80 * 1.1 / 0.9 * 0.99999]
    write_component_specbond(
        base,
        output,
        [_hit(aldehyde, amine_a, distances[0]), _hit(aldehyde, amine_b, distances[1])],
    )

    references = _reference_lengths(output, ("L4Y", "CE", "L5Y", "NZ"))
    assert references
    for distance in distances:
        assert any(0.9 * ref <= distance / 10.0 <= 1.1 * ref for ref in references)

    records, _ = read_specbond_records(output)
    matching = [
        record for record in records
        if (record[0], record[1], record[3], record[4]) == ("L4Y", "CE", "L5Y", "NZ")
    ]
    assert matching
    assert all(len(record[6].split(".")[1]) == 6 for record in matching)


def test_extra_crosslink_bond_fails_verification(tmp_path):
    gro = tmp_path / "component.gro"
    gro.write_text(
        "test\n4\n"
        + _gro_atom(1, "L4Y", "CE", 1, 0.000, 0.000, 0.000)
        + _gro_atom(2, "L5Y", "NZ", 2, 0.181, 0.000, 0.000)
        + _gro_atom(3, "L4Y", "CE", 3, 0.000, 1.000, 0.000)
        + _gro_atom(4, "L5Y", "NZ", 4, 0.181, 1.000, 0.000)
        + "1.0 1.0 1.0\n"
    )
    itp = tmp_path / "component.itp"
    itp.write_text(
        "[ atoms ]\n"
        "1 CT 1 L4Y CE 1 0.0 12.01\n"
        "2 N3 2 L5Y NZ 2 0.0 14.01\n"
        "3 CT 3 L4Y CE 3 0.0 12.01\n"
        "4 N3 4 L5Y NZ 4 0.0 14.01\n"
        "[ bonds ]\n1 2 1\n3 4 1\n"
    )
    assert count_crosslink_bonds(itp) == 2
    hit = _hit(
        _marker(0, "L4Y", "CE", 947, "A", (0.0, 0.0, 0.0)),
        _marker(1, "L5Y", "NZ", 9, "C", (1.81, 0.0, 0.0)),
        1.81,
    )
    with pytest.raises(RuntimeError, match="crosslink-like bonds"):
        verify_component_crosslinks(gro, itp, [hit])


def test_peptide_bond_is_not_counted_as_crosslink(tmp_path):
    itp = tmp_path / "peptide.itp"
    itp.write_text(
        "[ atoms ]\n"
        "1 CT 1 GLY C 1 0.0 12.01\n"
        "2 N3 2 ALA N 2 0.0 14.01\n"
        "[ bonds ]\n1 2 1\n"
    )
    assert count_crosslink_bonds(itp) == 0


def test_specbond_residues_from_crosslinks_csv_are_recognised():
    csv_path = (
        Path(__file__).parents[1] / "src/colbuilder/data/sequence/crosslinks.csv"
    )
    residues = set()
    with csv_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            for column in ("R1", "R2", "R3"):
                name = (row[column] or "").strip()
                if name and name != "NONE":
                    residues.add(name)
    assert {"DTY", "DTF"} <= residues

    bundled = (
        Path(__file__).parents[1]
        / "src/colbuilder/data/topology/amber99sb-star-ildnp.ff/specbond.dat"
    )
    records, _ = read_specbond_records(bundled)
    checked = []
    for record in records:
        if record[0] not in residues or record[3] not in residues:
            continue
        checked.append(record)
        assert record[1] in XLINK_ATOMS[record[0]]
        assert record[4] in XLINK_ATOMS[record[3]]
    assert len(checked) > 10
    assert any({record[0], record[3]} == {"DTY", "DTF"} for record in checked)


@pytest.mark.parametrize(
    "first, second, expected",
    [
        ("L4Y", "L5Y", True),
        ("L5Y", "L4Y", True),
        ("L4X", "L5X", True),
        ("LY4", "LY5", True),
        ("LX4", "LX5", True),
        ("LGX", "AGS", True),
        ("LPS", "APD", True),
        ("LZD", "LZS", True),
        ("DTY", "DTF", True),
        ("L4Y", "L5X", False),
        ("L4X", "L5Y", False),
        ("LY4", "L5Y", False),
        ("LGX", "APD", False),
        ("LPS", "AGS", False),
        ("L4Y", "L4Y", False),
    ],
)
def test_divalent_crosslinks_bond_only_matched_partners(first, second, expected):
    amber = Amber()
    first_link = Crosslink("1", first, "A", [0.0, 0.0, 0.0], "D", 0, "CE")
    second_link = Crosslink("1", second, "B", [1.8, 0.0, 0.0], "D", 1, "NZ")
    assert amber._are_compatible_crosslinks(first_link, second_link) is expected


def test_specbond_martini_allowed_without_topology():
    config = ColbuilderConfig(
        species="homo_sapiens",
        force_field="martini3",
        amber_topology_mode="specbond",
    )
    assert config.amber_topology_mode == "specbond"
