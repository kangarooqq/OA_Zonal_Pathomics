
from pathlib import Path
import ast


def test_zonal_depth_script_exists_and_parses():
    """Check that the zonal-depth module exists and has valid Python syntax."""
    script = Path(__file__).resolve().parents[1] / "simtk" / "zonal_depth.py"

    assert script.is_file(), "Missing simtk/zonal_depth.py"

    source = script.read_text(encoding="utf-8")
    ast.parse(source, filename=str(script))
