from pathlib import Path

root = Path(__file__).resolve().parents[1] / "apps" / "web"
for name in ("index.html", "app.js", "styles.css"):
    path = root / name
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"Missing or empty web asset: {path}")
html = (root / "index.html").read_text()
js = (root / "app.js").read_text()
for required in ("view-overview", "view-matrix", "view-strains", "view-genome", "view-families", "view-methods"):
    if required not in html:
        raise SystemExit(f"Missing view {required}")
if "function renderMatrix" not in js or "function drawDonut" not in js:
    raise SystemExit("Expected interactive visualization entry points are missing")
print("Frontend structure check passed.")
