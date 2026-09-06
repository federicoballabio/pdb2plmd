from __future__ import annotations
import csv
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src._version import VERSION

CONVERTER = ROOT / "pdb2plmd.py"
BASE = Path(__file__).resolve().parent
FIXTURES = BASE / "fixtures"
EXPECTED = BASE / "expected"
OUTPUT = BASE / "output"
CASE_TABLE = next(BASE.glob("*.csv"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def readiness_from_log(path: Path) -> str:
    if not path.exists():
        return ""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ONEBEAD_SAXS_READY "):
            return line.split(None, 1)[1].strip()
    return ""


def compact(text: str, limit: int = 6) -> str:
    lines = [line for line in text.strip().splitlines() if line.strip()]
    if not lines:
        return "<empty>"
    if len(lines) <= limit:
        return "\n".join(lines)
    return "\n".join(lines[:limit] + [f"... ({len(lines) - limit} more line(s); see generated log)"])


def load_cases():
    with CASE_TABLE.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def command_text(pdb_id: str, verbose: bool) -> str:
    parts = ["python3", "pdb2plmd.py", "-i", f"examples/fixtures/{pdb_id}.pdb", "-o", f"examples/output/{pdb_id}_out.pdb", "--drop-solvent"]
    if verbose:
        parts.append("-v")
    return " ".join(parts)


def main() -> int:
    cases = load_cases()
    if OUTPUT.exists():
        shutil.rmtree(OUTPUT)
    OUTPUT.mkdir()
    passed = 0
    for index, case in enumerate(cases, 1):
        pdb_id = case["pdb_id"]
        verbose = case["verbose"] == "1"
        expected_rc = int(case["expected_rc"])
        inp = FIXTURES / f"{pdb_id}.pdb"
        out = OUTPUT / f"{pdb_id}_out.pdb"
        log = out.with_suffix(".log")
        expected_pdb = EXPECTED / f"{pdb_id}_out.pdb"
        fixture_ok = inp.exists() and sha256(inp) == case["fixture_sha256"]
        cmd = [sys.executable, str(CONVERTER), "-i", str(inp), "-o", str(out), "--drop-solvent"]
        if verbose:
            cmd.append("-v")
        proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        readiness = readiness_from_log(log)
        ok = fixture_ok and proc.returncode == expected_rc
        if expected_rc == 0:
            ok = ok and out.exists()
            if verbose:
                ok = ok and log.exists() and proc.stdout.startswith(f"pdb2plmd {VERSION}: OK - wrote ")
                ok = ok and readiness == case["expected_readiness"]
            else:
                ok = ok and not log.exists() and proc.stdout == f"pdb2plmd {VERSION}: PDB converted.\n"
            ok = ok and expected_pdb.exists() and out.read_bytes() == expected_pdb.read_bytes()
            ok = ok and sha256(out) == case["expected_pdb_sha256"]
        else:
            ok = ok and not out.exists() and log.exists()
            marker = case["error_marker"]
            log_text = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
            ok = ok and marker in log_text
            if verbose:
                ok = ok and proc.stderr.startswith(f"pdb2plmd {VERSION}: ERROR - ")
            else:
                ok = ok and proc.stderr == f"pdb2plmd {VERSION}: PDB conversion failed.\n"
        expected_text = "PDB converted" if expected_rc == 0 else "PDB conversion failed"
        if case["expected_readiness"]:
            expected_text += f"; ONEBEAD_SAXS_READY={case['expected_readiness']}"
        print(f"[{index}/{len(cases)}] {pdb_id}: {case['label']}")
        print("command:", command_text(pdb_id, verbose))
        print("expected:", expected_text)
        print("fixture SHA256:", "PASS" if fixture_ok else "FAIL")
        print("observed return code:", proc.returncode)
        if readiness:
            print("observed ONEBEAD_SAXS_READY:", readiness)
        print("observed stdout:", compact(proc.stdout))
        print("observed stderr:", compact(proc.stderr))
        if expected_rc == 0:
            print("verified output PDB:", "PASS" if out.exists() and expected_pdb.exists() and out.read_bytes() == expected_pdb.read_bytes() else "FAIL")
        print("verified:", "PASS" if ok else "FAIL")
        print()
        passed += int(ok)
    print(f"RESULT: {passed}/{len(cases)} PASS")
    return 0 if passed == len(cases) else 1


if __name__ == "__main__":
    raise SystemExit(main())
