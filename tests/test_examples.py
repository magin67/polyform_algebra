"""Keep every documented example runnable and safe to import."""
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = sorted(path.stem for path in (ROOT / "examples").glob("*.py")
                  if path.stem != "__init__")


@pytest.mark.parametrize("name", EXAMPLES)
def test_example_runs_as_module(name):
    """Catch removed APIs and invalid input syntax in the documented module command."""
    result = subprocess.run(
        [sys.executable, "-B", "-m", f"examples.{name}"],
        cwd=ROOT, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.strip(), "The example should display its calculation."


@pytest.mark.parametrize("name", EXAMPLES)
def test_example_import_has_no_output(name):
    """Importing a demonstration must not start calculations or print results."""
    result = subprocess.run(
        [sys.executable, "-B", "-c", f"import examples.{name}"],
        cwd=ROOT, capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""
