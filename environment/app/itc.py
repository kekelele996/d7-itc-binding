"""批次复核的快速换算脚本。

给出一个批次目录，按首几针的热量估结合焓变，按热量掉下来的针位估化学计量比，
再按半饱和处的游离配体浓度估结合常数。结果从标准输出打一行 JSON。
"""

import json
import math
import sys

import itc_io


def _plateau_heat(heats):
    """首几针的平均热量，作为饱和区的热量水平。"""
    head = heats[:3]
    return sum(head) / len(head)


def estimate(batch):
    blank_runs, sample_runs = itc_io.split_roles(batch)
    if not sample_runs:
        raise SystemExit("批次里没有样品滴定")
    run = sample_runs[0]

    v0 = run["cell_volume_uL"]
    cell = run["cell_concentration_M"]
    syringe = run["syringe_concentration_M"]

    heats = [item["heat_uJ"] for item in run["injections"]]
    volumes = [item["volume_uL"] for item in run["injections"]]

    plateau = _plateau_heat(heats)
    # 首针注进去的配体全部结合时的热量，用饱和区的热量水平反推焓变。
    first_ligand = syringe * volumes[0] / v0
    if abs(plateau) < 1e-9 or first_ligand <= 0.0:
        raise SystemExit("首针热量太低，无法换算")
    dh_j_per_mol = plateau / (v0 * 1e-6 * first_ligand * 1e6)

    # 热量掉到饱和区两成的位置当作等当点。
    threshold = 0.2 * abs(plateau)
    equivalence = len(heats)
    for idx, heat in enumerate(heats):
        if abs(heat) <= threshold:
            equivalence = idx
            break

    ligand_total = 0.0
    for idx in range(min(equivalence, len(volumes))):
        ligand_total += syringe * volumes[idx] / v0
    n = ligand_total / cell if cell > 0 else 0.0

    # 半饱和处的游离配体浓度取等当点之前一半的累计配体。
    half = ligand_total / 2.0
    free_ligand = max(half, 1e-12)
    k = 1.0 / free_ligand if free_ligand > 0 else 0.0

    return {"batch_id": batch["batch_id"],
            "K_M-1": k,
            "dH_kJ_per_mol": dh_j_per_mol / 1000.0,
            "n": n}


def main(argv):
    if len(argv) != 2:
        return 2
    batch = itc_io.load_batch(argv[1])
    result = estimate(batch)
    result["K_M-1"] = float(f"{result['K_M-1']:.3g}")
    result["dH_kJ_per_mol"] = round(result["dH_kJ_per_mol"], 2)
    result["n"] = round(result["n"], 3)
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
