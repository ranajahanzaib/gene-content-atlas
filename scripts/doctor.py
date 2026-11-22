import importlib.util
import shutil
import subprocess
import sys


def version(label: str, executable: str, args: list[str]) -> None:
    path = shutil.which(executable)
    if not path:
        print(f"{label}: not installed")
        return
    try:
        result = subprocess.run([path, *args], text=True, capture_output=True, timeout=15, check=False)
    except subprocess.TimeoutExpired:
        print(f"{label}: installed at {path}; version check timed out")
        return
    output = (result.stdout or result.stderr).strip().splitlines()
    print(f"{label}: {output[0] if output else f'installed at {path}, version unavailable'}")


print(f"Python: {sys.version.split()[0]}")
version("Node.js (optional)", "node", ["--version"])
version("Docker (optional)", "docker", ["--version"])
version("NCBI Datasets CLI", "datasets", ["--version"])
version("Panaroo", "panaroo", ["--version"])
version("MMseqs2", "mmseqs", ["version"])
version("MAFFT", "mafft", ["--version"])
version("IQ-TREE", "iqtree2", ["--version"])
version("Redis", "redis-server", ["--version"])
print("API dependencies: Python standard library only")
print(f"Biopython: {'installed' if importlib.util.find_spec('Bio') else 'not installed (not required by current parsers)'}")
print("Database: not configured (current artifact store is filesystem-backed)")
for directory in ("data", "data/real", "data/generated"):
    print(f"Storage directory {directory}: {'ready' if __import__('pathlib').Path(directory).exists() else 'will be created when needed'}")
print("NCBI downloads and Panaroo execution require their respective executables and network access.")
