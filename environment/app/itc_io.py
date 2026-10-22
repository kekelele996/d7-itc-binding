"""批次目录的读取与汇总（仪器导出格式的适配层）。"""

import json
from pathlib import Path


def _read_injections(path):
    rows = []
    text = path.read_text(encoding="utf-8")
    for line in text.splitlines()[1:]:
        line = line.strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",")]
        rows.append({
            "injection": int(parts[0]),
            "volume_uL": float(parts[1]),
            "heat_uJ": float(parts[2]),
        })
    rows.sort(key=lambda item: item["injection"])
    return rows


def load_batch(batch_dir):
    """读取批次目录，返回 {"batch_id", "sample", "runs": [...]}。"""
    root = Path(batch_dir)
    meta = json.loads((root / "batch.json").read_text(encoding="utf-8"))
    runs = []
    runs_dir = root / "runs"
    for run_dir in sorted(p for p in runs_dir.iterdir() if p.is_dir()):
        run_meta = json.loads((run_dir / "meta.json").read_text(encoding="utf-8"))
        injections = _read_injections(run_dir / "injections.csv")
        runs.append({
            "run_id": run_meta.get("run_id", run_dir.name),
            "role": run_meta.get("role", "sample"),
            "cell_volume_uL": float(run_meta["cell_volume_uL"]),
            "cell_concentration_M": float(run_meta["cell_concentration_M"]),
            "syringe_concentration_M": float(run_meta["syringe_concentration_M"]),
            "temperature_C": float(run_meta.get("temperature_C", 25.0)),
            "injections": injections,
        })
    return {"batch_id": meta.get("batch_id", root.name),
            "sample": meta.get("sample", ""),
            "note": meta.get("note", ""),
            "runs": runs}


def split_roles(batch):
    blank = [r for r in batch["runs"] if r["role"] == "blank"]
    samples = [r for r in batch["runs"] if r["role"] != "blank"]
    return blank, samples
