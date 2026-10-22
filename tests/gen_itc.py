"""隐藏批次的生成器：从给定的结合参数正演出滴定记录。

平衡用二分法解质量守恒方程，和参考实现里的二次方程闭式解法是两条独立路径。
生成参数就是期望值，不来自参考解的输出。
"""

import csv
import json
import math
import random
from pathlib import Path

NOISE_UI = 0.2
CELL_DH = 8e-6
DILUTION_DH = 2.5


def _free_ligand(cell_conc, total_ligand, k, n):
    """二分求游离配体浓度；结合量随游离浓度单调增。"""
    if cell_conc <= 0.0:
        return total_ligand
    lo, hi = 0.0, total_ligand
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        bound = n * cell_conc * k * mid / (1.0 + k * mid)
        if mid + bound > total_ligand:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def simulate_heats(system, protocol, dilution_per_ul, noise, rng):
    """返回 [(volume_uL, heat_uJ)]，含池体积置换和混合热。"""
    v0 = protocol["cell_volume_uL"]
    cell = protocol["cell_concentration_M"]
    syringe = protocol["syringe_concentration_M"]
    total_ligand = 0.0
    prev_bound = 0.0
    rows = []
    for volume in protocol["volumes_uL"]:
        frac = volume / v0
        cell -= cell * frac
        total_ligand = total_ligand * (1.0 - frac) + syringe * frac
        free = _free_ligand(cell, total_ligand, system["K"], system["n"])
        bound = system["n"] * cell * system["K"] * free / (1.0 + system["K"] * free)
        heat = v0 * 1e-6 * system["dH_J"] * (bound - prev_bound * (1.0 - frac)) * 1e6
        prev_bound = bound
        rows.append((volume, heat + dilution_per_ul * volume + rng.gauss(0.0, noise)))
    return rows


def _write_run(root, run_id, role, protocol, rows):
    run_dir = root / "runs" / run_id
    run_dir.mkdir(parents=True)
    meta = {"run_id": run_id, "role": role,
            "cell_volume_uL": protocol["cell_volume_uL"],
            "cell_concentration_M": protocol["cell_concentration_M"],
            "syringe_concentration_M": protocol["syringe_concentration_M"],
            "temperature_C": 25.0}
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    with (run_dir / "injections.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["injection", "volume_uL", "heat_uJ"])
        for idx, (volume, heat) in enumerate(rows, start=1):
            writer.writerow([idx, f"{volume:.2f}", f"{heat:.4f}"])


# ---------- 三组批次 ----------

def _dh_batch(system, seed):
    rng = random.Random(seed)
    volumes = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4]
    cell = CELL_DH
    syringe = 2.5 * cell * 800.0 / sum(volumes)
    sample = {"cell_volume_uL": 800.0, "cell_concentration_M": cell,
              "syringe_concentration_M": syringe, "volumes_uL": volumes}
    blank = dict(sample, cell_concentration_M=0.0)
    return ([("S1", "sample", sample, simulate_heats(system, sample, DILUTION_DH, NOISE_UI, rng))],
            [("B1", "blank", blank,
              [(v, DILUTION_DH * v + rng.gauss(0.0, NOISE_UI)) for v in volumes])])


def _n_batch(system, seed):
    rng = random.Random(seed)
    volumes = [2.0] * 25
    cell = 6.0e-5
    syringe = 2.5 * cell * 200.0 / sum(volumes)
    sample = {"cell_volume_uL": 200.0, "cell_concentration_M": cell,
              "syringe_concentration_M": syringe, "volumes_uL": volumes}
    return ([("S1", "sample", sample, simulate_heats(system, sample, 0.0, NOISE_UI, rng))], [])


def _bad_batch(system, seed):
    rng = random.Random(seed)
    volumes = [1.0] * 20
    cell = 7.0 / (system["K"] * system["n"])
    syringe = 2.5 * cell * 1000.0 / sum(volumes)
    sample = {"cell_volume_uL": 1000.0, "cell_concentration_M": cell,
              "syringe_concentration_M": syringe, "volumes_uL": volumes}
    rows = simulate_heats(system, sample, 0.0, 0.15, rng)
    rows[7] = (rows[7][0], rows[7][1] + 20.0)   # 进样气泡
    rows[15] = (rows[15][0], rows[15][1] + 9.0)  # 溅射
    return ([("S1", "sample", sample, rows)], [])


BUILDERS = {"dh": _dh_batch, "n": _n_batch, "bad": _bad_batch}

SYSTEMS = {
    "dh": [{"K": 1.5e6, "dH_J": -45000.0, "n": 1.00, "seed": 41201},
           {"K": 8.0e5, "dH_J": -38000.0, "n": 1.30, "seed": 41202},
           {"K": 2.2e6, "dH_J": -52000.0, "n": 0.90, "seed": 41203}],
    "n": [{"K": 3.0e5, "dH_J": -30000.0, "n": 1.20, "seed": 41301},
          {"K": 5.0e5, "dH_J": -42000.0, "n": 0.95, "seed": 41302},
          {"K": 1.8e5, "dH_J": -35000.0, "n": 1.40, "seed": 41303}],
    "bad": [{"K": 1.0e6, "dH_J": -50000.0, "n": 1.00, "seed": 41401},
            {"K": 4.0e5, "dH_J": -44000.0, "n": 1.10, "seed": 41402},
            {"K": 2.0e6, "dH_J": -56000.0, "n": 0.90, "seed": 41403}],
}


def expected(kind, index):
    system = SYSTEMS[kind][index]
    return {"K_M-1": system["K"], "dH_kJ_per_mol": system["dH_J"] / 1000.0,
            "n": system["n"]}


def build_batch(kind, index, dest):
    """在 dest 下写一个批次目录，返回 (batch_id, 期望值)。"""
    system = SYSTEMS[kind][index]
    sample_runs, blank_runs = BUILDERS[kind](system, system["seed"])
    root = Path(dest) / f"B-{kind}-{index + 1}"
    root.mkdir(parents=True, exist_ok=True)
    batch_id = f"B-{kind}-{index + 1}"
    (root / "batch.json").write_text(json.dumps(
        {"batch_id": batch_id, "sample": "XRD-CER-071",
         "note": "复核批次"}, indent=2), encoding="utf-8")
    for run_id, role, protocol, rows in sample_runs + blank_runs:
        _write_run(root, run_id, role, protocol, rows)
    return batch_id, expected(kind, index)
