"""End-to-end smoke test on the synthetic dataset (needs the external tools).

Point PHILOMENOS_TOOLS_CONFIG at a YAML file with a `tools:` section, e.g.
    PHILOMENOS_TOOLS_CONFIG=config/tools.local.yaml pytest -m smoke
Without it the test is skipped. It is equivalent to `philomenos smoke`.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest


@pytest.mark.smoke
def test_smoke_run(tmp_path: Path) -> None:
    tools = os.environ.get("PHILOMENOS_TOOLS_CONFIG")
    if not tools:
        pytest.skip("set PHILOMENOS_TOOLS_CONFIG to a config with a tools: section")
    from philomenos.smoke import run_smoke
    assert run_smoke(tmp_path, tools_config=Path(tools), keep=False) == 0
