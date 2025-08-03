from pathlib import Path
root   = Path(__file__).resolve().parents[1]
data   = root / "data"
tidy   = data / "tidy"
splits = data / "splits"
models = root / "models"
