"""批次复核：从热量-注射记录反演结合常数、结合焓变和化学计量比。

正演：池内物质随每次注射按体积分数被带走，配体按注射量补入，结合平衡按
1:1 多位点解二次方程；每针热量正比于「本针结合量 − 上一针残留结合量」。
反演：把 dH 做解析消元，只对 (K, n) 做网格加局部细化，再用迭代重加权压掉
异常注射点。
"""

import json
import math
import sys

import itc_io

SUBTRACT_BLANK = True
MODEL_DISPLACEMENT = True
ROBUST_WEIGHTS = True

DH_MIN = 1.0e2
DH_MAX = 2.0e6


def bound_ligand(cell_conc, total_ligand, k, n):
    """平衡时已结合配体的浓度（与 cell_conc 同单位）。"""
    if cell_conc <= 0.0:
        return 0.0
    a = 1.0 + k * n * cell_conc - k * total_ligand
    free = (-a + math.sqrt(a * a + 4.0 * k * total_ligand)) / (2.0 * k)
    return n * cell_conc * k * free / (1.0 + k * free)


def unit_heats(k, n, run, displacement):
    """每针的单位焓热量：乘上 dH（J/mol）就是 µJ。"""
    v0 = run["cell_volume_uL"]
    cell = run["cell_concentration_M"]
    syringe = run["syringe_concentration_M"]
    total_ligand = 0.0
    prev = 0.0
    out = []
    for item in run["injections"]:
        frac = item["volume_uL"] / v0
        if displacement:
            cell *= (1.0 - frac)
            total_ligand = total_ligand * (1.0 - frac) + syringe * frac
            bound = bound_ligand(cell, total_ligand, k, n)
            out.append(v0 * 1e-6 * (bound - prev * (1.0 - frac)) * 1e6)
        else:
            total_ligand += syringe * frac
            bound = bound_ligand(cell, total_ligand, k, n)
            out.append(v0 * 1e-6 * (bound - prev) * 1e6)
        prev = bound
    return out


def dilution_slope(blank_runs):
    """从空白滴定回归每微升的混合热（µJ/µL），没有空白就取 0。"""
    if not blank_runs:
        return 0.0
    num = 0.0
    den = 0.0
    for run in blank_runs:
        for item in run["injections"]:
            v = item["volume_uL"]
            num += v * item["heat_uJ"]
            den += v * v
    return num / den if den > 0 else 0.0


def build_datasets(batch, slope):
    """把每条样品滴定的热量换成扣掉混合热之后的观测值。"""
    datasets = []
    for run in batch["runs"]:
        if run["role"] == "blank":
            continue
        heats = [item["heat_uJ"] for item in run["injections"]]
        if SUBTRACT_BLANK:
            heats = [h - slope * item["volume_uL"]
                     for h, item in zip(heats, run["injections"])]
        datasets.append((run, heats))
    return datasets


def _weighted_cost(k, n, datasets, weights):
    uu = uq = 0.0
    units = []
    pos = 0
    for run, heats in datasets:
        unit = unit_heats(k, n, run, MODEL_DISPLACEMENT)
        units.append(unit)
        for ui, qi in zip(unit, heats):
            w = weights[pos]
            pos += 1
            uu += w * ui * ui
            uq += w * ui * qi
    if uu <= 0.0:
        return None
    dh = uq / uu
    if not (DH_MIN <= abs(dh) <= DH_MAX):
        return None
    cost = 0.0
    pos = 0
    for (run, heats), unit in zip(datasets, units):
        for ui, qi in zip(unit, heats):
            r = qi - dh * ui
            cost += weights[pos] * r * r
            pos += 1
    return cost, dh


def _search(datasets, weights):
    """两维粗网格 + 三级局部细化。"""
    best = None
    for i in range(41):
        k = 10.0 ** (2.0 + 8.0 * i / 40.0)
        for j in range(39):
            n = 0.3 + 1.9 * j / 38.0
            got = _weighted_cost(k, n, datasets, weights)
            if got and (best is None or got[0] < best[0]):
                best = (got[0], k, got[1], n)
    if best is None:
        return None
    for span_k, span_n, steps in ((0.15, 0.04, 16), (0.015, 0.004, 16), (0.0015, 0.0004, 10)):
        logk0 = math.log10(best[1])
        n0 = best[3]
        for di in range(-steps, steps + 1):
            k = 10.0 ** (logk0 + span_k * di / steps)
            for dj in range(-steps, steps + 1):
                n = n0 + span_n * dj / steps
                if n <= 0.01:
                    continue
                got = _weighted_cost(k, n, datasets, weights)
                if got and got[0] < best[0]:
                    best = (got[0], k, got[1], n)
    return {"K": best[1], "dH": best[2], "n": best[3]}


def residuals(fit, datasets):
    out = []
    for run, heats in datasets:
        unit = unit_heats(fit["K"], fit["n"], run, MODEL_DISPLACEMENT)
        out.extend(qi - fit["dH"] * ui for ui, qi in zip(unit, heats))
    return out


def fit_parameters(datasets):
    """先做一次等权拟合；按残差的稳健尺度做迭代重加权。"""
    npts = sum(len(heats) for _, heats in datasets)
    weights = [1.0] * npts
    fit = None
    rounds = 6 if ROBUST_WEIGHTS else 1
    for _ in range(rounds):
        fit = _search(datasets, weights)
        if fit is None:
            return None
        if not ROBUST_WEIGHTS:
            break
        res = residuals(fit, datasets)
        scale = sorted(abs(r) for r in res)[len(res) // 2]
        sigma = 1.4826 * scale if scale > 1e-12 else 1e-9
        weights = [1.0 / (1.0 + (abs(r) / (2.5 * sigma)) ** 2) for r in res]
    return fit


def analyse(batch):
    slope = dilution_slope([r for r in batch["runs"] if r["role"] == "blank"])
    datasets = build_datasets(batch, slope)
    if not datasets:
        raise SystemExit("批次里没有样品滴定")
    fit = fit_parameters(datasets)
    if fit is None:
        raise SystemExit("拟合没有收敛")
    return {"batch_id": batch["batch_id"],
            "K_M-1": float(f"{fit['K']:.3g}"),
            "dH_kJ_per_mol": round(fit["dH"] / 1000.0, 2),
            "n": round(fit["n"], 3)}


def main(argv):
    if len(argv) != 2:
        return 2
    batch = itc_io.load_batch(argv[1])
    result = analyse(batch)
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
