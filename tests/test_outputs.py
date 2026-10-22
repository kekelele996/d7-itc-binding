"""结合参数复检的验收测试。

隐藏批次由 gen_itc 现造：期望值就是生成时用的结合参数，不读参考解的输出。
每个 test_ 的 docstring 用「」引用题面里的一句话；测试体里不写 if/try，
分支和异常处理放在模块级辅助函数里。被测程序一律以 runner 身份、干净环境变量、
独立进程组运行，超时整组结束。
"""

import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, "/tests")
import gen_itc  # noqa: E402

ENTRY = "/app/itc.py"
PYTHON = "/usr/local/bin/python3"
CANDIDATE_USER = "runner"
CANDIDATE_ENV = {
    "PATH": "/usr/local/bin:/usr/bin:/bin",
    "HOME": "/home/runner",
    "LANG": "C.UTF-8",
    "TZ": "UTC",
}
CALL_TIMEOUT_SEC = 60


def run_candidate(*args, timeout=CALL_TIMEOUT_SEC):
    """以 runner 身份运行候选入口；超时结束整个进程组并抛出 TimeoutError。"""
    as_root = os.getuid() == 0
    proc = subprocess.Popen(
        [PYTHON, ENTRY, *[str(a) for a in args]],
        cwd="/app",
        env=CANDIDATE_ENV,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
        user=CANDIDATE_USER if as_root else None,
        group=CANDIDATE_USER if as_root else None,
        extra_groups=[] if as_root else None,
    )
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.communicate()
        raise TimeoutError(f"{ENTRY} 超过 {timeout}s")
    return subprocess.CompletedProcess(proc.args, proc.returncode, out, err)


def _open_permissions(root):
    """批次目录要给 runner 读：目录 0755、文件 0644。"""
    root.chmod(0o755)
    for path in root.rglob("*"):
        path.chmod(0o755 if path.is_dir() else 0o644)


def analyse(kind, index):
    """造一个隐藏批次，跑候选，返回 (输出 JSON, 期望参数)。"""
    workdir = Path(tempfile.mkdtemp(prefix="itc-"))
    batch_id, expected = gen_itc.build_batch(kind, index, workdir)
    root = workdir / batch_id
    _open_permissions(workdir)
    _open_permissions(root)
    proc = run_candidate(root)
    result = {"returncode": proc.returncode, "stdout": proc.stdout}
    shutil.rmtree(workdir, ignore_errors=True)
    assert result["returncode"] == 0, f"退出码 {result['returncode']}"
    return json.loads(result["stdout"]), expected


def relative_error(got, expected, key):
    return abs(got[key] - expected[key]) / abs(expected[key])


# ---------- gate：入口与输出契约，不计分 ----------

def test_prints_json_object():
    """「标准输出一个 JSON 对象，只有四个字段」"""
    got, _ = analyse("n", 0)
    assert set(got) == {"batch_id", "K_M-1", "dH_kJ_per_mol", "n"}


def test_wrong_argument_count_exits_2():
    """「命令行参数不对时直接以退出码 2 结束」"""
    proc = run_candidate()
    assert proc.returncode == 2
    assert proc.stdout.strip() == ""


# ---------- 难点 1：混合热与空白滴定 ----------

def test_enthalpy_batch1():
    """「dH 的相对误差不超过 2%」"""
    got, expected = analyse("dh", 0)
    assert relative_error(got, expected, "dH_kJ_per_mol") <= 0.02


def test_enthalpy_batch2():
    """「dH 的相对误差不超过 2%」"""
    got, expected = analyse("dh", 1)
    assert relative_error(got, expected, "dH_kJ_per_mol") <= 0.02


def test_enthalpy_batch3():
    """「dH 的相对误差不超过 2%」"""
    got, expected = analyse("dh", 2)
    assert relative_error(got, expected, "dH_kJ_per_mol") <= 0.02


# ---------- 难点 2：池体积置换 ----------

def test_stoichiometry_batch1():
    """「n 的相对误差不超过 4%」"""
    got, expected = analyse("n", 0)
    assert relative_error(got, expected, "n") <= 0.04


def test_stoichiometry_batch2():
    """「n 的相对误差不超过 4%」"""
    got, expected = analyse("n", 1)
    assert relative_error(got, expected, "n") <= 0.04


def test_stoichiometry_batch3():
    """「n 的相对误差不超过 4%」"""
    got, expected = analyse("n", 2)
    assert relative_error(got, expected, "n") <= 0.04


# ---------- 难点 3：异常注射点的稳健性 ----------

def test_abnormal_injections_batch1():
    """「偶发的进样异常也照原样留着」"""
    got, expected = analyse("bad", 0)
    assert relative_error(got, expected, "K_M-1") <= 0.10
    assert relative_error(got, expected, "dH_kJ_per_mol") <= 0.03
    assert relative_error(got, expected, "n") <= 0.04


def test_abnormal_injections_batch2():
    """「偶发的进样异常也照原样留着」"""
    got, expected = analyse("bad", 1)
    assert relative_error(got, expected, "K_M-1") <= 0.10
    assert relative_error(got, expected, "dH_kJ_per_mol") <= 0.03
    assert relative_error(got, expected, "n") <= 0.04


def test_abnormal_injections_batch3():
    """「偶发的进样异常也照原样留着」"""
    got, expected = analyse("bad", 2)
    assert relative_error(got, expected, "K_M-1") <= 0.10
    assert relative_error(got, expected, "dH_kJ_per_mol") <= 0.03
    assert relative_error(got, expected, "n") <= 0.04
