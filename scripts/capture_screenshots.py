"""Automated screenshot capture for Gene Content Atlas."""
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCREENSHOTS_DIR = ROOT / "docs" / "screenshots"
SCREENSHOTS_DIR.mkdir(parents=True, exist_ok=True)

CHROME_PATH = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
PORT = 8765
BASE_URL = f"http://127.0.0.1:{PORT}"

SHOTS = [
    # Light Mode
    ("overview.png", f"{BASE_URL}/#overview", 1440, 920),
    ("matrix.png", f"{BASE_URL}/#matrix", 1440, 920),
    ("strains.png", f"{BASE_URL}/#strains", 1440, 900),
    ("genome-landscape.png", f"{BASE_URL}/#genome", 1440, 900),
    ("families.png", f"{BASE_URL}/#families", 1440, 900),
    ("tree.png", f"{BASE_URL}/?strains=LP-01,LP-02,LP-03,LP-04#tree", 1440, 900),
    ("statistics.png", f"{BASE_URL}/#statistics", 1440, 980),
    ("compare.png", f"{BASE_URL}/?strains=LP-01,LP-02,LP-03#compare", 1440, 920),
    ("real-dataset-overview.png", f"{BASE_URL}/?dataset=real-lp-refseq-panaroo#overview", 1440, 920),
    ("real-dataset-matrix.png", f"{BASE_URL}/?dataset=real-lp-refseq-panaroo#matrix", 1440, 920),

    # Dark Mode
    ("overview-dark.png", f"{BASE_URL}/?theme=dark#overview", 1440, 920),
    ("matrix-dark.png", f"{BASE_URL}/?theme=dark#matrix", 1440, 920),
    ("tree-dark.png", f"{BASE_URL}/?theme=dark&strains=LP-01,LP-02,LP-03,LP-04#tree", 1440, 900),
    ("statistics-dark.png", f"{BASE_URL}/?theme=dark#statistics", 1440, 980),
    ("compare-dark.png", f"{BASE_URL}/?theme=dark&strains=LP-01,LP-02,LP-03#compare", 1440, 920),
    ("real-dataset-overview-dark.png", f"{BASE_URL}/?theme=dark&dataset=real-lp-refseq-panaroo#overview", 1440, 920),
]

def capture(only_dark=False):
    print(f"Capturing screenshots from {BASE_URL}...")
    shots = [s for s in SHOTS if ("-dark" in s[0]) == only_dark] if only_dark else SHOTS
    for filename, url, width, height in shots:
        out_path = SCREENSHOTS_DIR / filename
        cmd = [
            CHROME_PATH,
            "--headless=new",
            "--disable-gpu",
            f"--window-size={width},{height}",
            "--virtual-time-budget=3500",
            "--hide-scrollbars",
            f"--screenshot={out_path}",
            url,
        ]
        print(f"Capturing {filename} ({width}x{height}) from {url}...")
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            print(f"Error capturing {filename}: {res.stderr}")
        else:
            size_kb = out_path.stat().st_size / 1024 if out_path.exists() else 0
            print(f"Saved {filename} ({size_kb:.1f} KB)")

if __name__ == "__main__":
    only_dark = "--only-dark" in sys.argv
    capture(only_dark=only_dark)
