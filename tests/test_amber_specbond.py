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
from colbuilder.core.topology.amber_specbond import (
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
    itp.write_text("[ bonds ]\n1 2 1\n")
    hit = _hit(
        _marker(0, "L4Y", "CE", 947, "A", (0.0, 0.0, 0.0)),
        _marker(1, "L5Y", "NZ", 9, "C", (1.81, 0.0, 0.0)),
        1.81,
    )
    assert verify_component_crosslinks(gro, itp, [hit]) == 1

    itp.write_text("[ bonds ]\n")
    with pytest.raises(RuntimeError, match="did not create expected crosslink bond"):
        verify_component_crosslinks(gro, itp, [hit])
