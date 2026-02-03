#!/usr/bin/env python3
"""
Topology helper for mixed fibrils.

Usage:
  python scripts/topology_from_mix.py <config.yaml>
"""

from __future__ import annotations

import asyncio
import shutil
from pathlib import Path
from typing import Dict, List

from colbuilder.core.geometry.model import Model
from colbuilder.core.geometry.system import System
from colbuilder.core.topology.martini import build_martini3
from colbuilder.core.utils.config import load_yaml_config, validate_config
from colbuilder.core.utils.files import FileManager
from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.martinize_finder import find_and_install_custom_force_field


LOG = setup_logger(__name__)


def _parse_connect_file(connect_file: Path) -> Dict[float, Dict[str, List[float]]]:
    """
    Parse connect_from_colbuilder.txt into a mapping:
      model_id -> {"type": "A", "connect": [ids...]}
    """
    mapping: Dict[float, Dict[str, List[float]]] = {}
    with open(connect_file, "r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if ";" in line:
                left, right = line.split(";", 1)
                model_type = right.strip()
            else:
                left, model_type = line, "A"
            ids = [
                float(tok.replace(".caps.pdb", ""))
                for tok in left.split()
                if tok.strip()
            ]
            if not ids:
                continue
            for model_id in ids:
                mapping[model_id] = {"type": model_type, "connect": ids}
    return mapping


def _build_system_from_connect(connect_file: Path) -> System:
    data = _parse_connect_file(connect_file)
    system = System()
    for model_id, info in data.items():
        model = Model(id=model_id, transformation=[0.0, 0.0, 0.0])
        model.type = info["type"]
        model.add_connect(connect_id=model_id, connect=info["connect"])
        system.add_model(model)
    return system


def _copy_caps(src_root: Path, dest_root: Path) -> None:
    for cap_file in src_root.glob("**/*.caps.pdb"):
        rel_path = cap_file.relative_to(src_root)
        dest_dir = dest_root / rel_path.parent
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(cap_file, dest_dir / cap_file.name)


async def _run(config_path: Path) -> None:
    cfg_dict = load_yaml_config(config_path)
    cfg = validate_config(cfg_dict)

    if not cfg.force_field or cfg.force_field != "martini3":
        raise ValueError("force_field must be 'martini3' for this helper.")

    working_dir = Path(cfg.working_directory).resolve()
    mixing_dir = working_dir / ".tmp" / "mixing_crosslinks"
    connect_file = mixing_dir / "connect_from_colbuilder.txt"
    if not connect_file.exists():
        raise FileNotFoundError(f"connect file not found: {connect_file}")

    system = _build_system_from_connect(connect_file)

    file_manager = FileManager(cfg)
    topology_dir = file_manager.get_temp_dir("topology_gen")
    topology_dir.mkdir(parents=True, exist_ok=True)

    if not mixing_dir.exists():
        raise FileNotFoundError(f"mixing directory not found: {mixing_dir}")

    _copy_caps(mixing_dir, topology_dir)

    find_and_install_custom_force_field(cfg.FORCE_FIELD_DIR)
    await build_martini3(system=system, config=cfg, file_manager=file_manager)


def main() -> None:
    import sys

    if len(sys.argv) != 2:
        print("Usage: python scripts/topology_from_mix.py <config.yaml>")
        sys.exit(1)

    config_path = Path(sys.argv[1]).resolve()
    asyncio.run(_run(config_path))


if __name__ == "__main__":
    main()
