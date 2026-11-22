from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "services" / "api"))
from app import DATA, demo_data  # noqa: E402

DATA.parent.mkdir(parents=True, exist_ok=True)
DATA.write_text(__import__("json").dumps(demo_data(), indent=2) + "\n")
print(f"Wrote deterministic synthetic fixture: {DATA}")
