#!/usr/bin/env python3
"""
00_verify_environment_lock_(Ryan).py

Check this machine against ENVIRONMENT_LOCK.md *before* running anything expensive.

Why this exists as a script rather than a checklist: the HADDOCK3 venv is addressed by an
absolute path baked in at venv-creation time, in the shebang of every console script. Move the
project folder and 05/06/08/09 stop launching -- and the failure surfaces as "bad interpreter",
deep inside a subprocess, looking like a pipeline bug rather than an environment one. This is
the check that names the broken layer in two seconds.

It also refuses to let an architecture mismatch pass silently: PyMOL.app here is x86_64 on an
arm64 host.

Stdlib only, so it runs before anything has been repaired or installed.

Run:
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3 \
        "00_verify_environment_lock_(Ryan).py"

Exit: 0 = every check passed. 1 = at least one FAIL (WARN alone does not fail).
"""

from __future__ import annotations

import json
import platform
import plistlib
import subprocess
import sys
from pathlib import Path

# --------------------------------------------------------------------------- paths under test

PROJECT = Path(__file__).resolve().parent
if not (PROJECT / "PROTOCOL_LOCK.md").exists():          # tolerate being called from elsewhere
    PROJECT = Path.cwd()

SYSTEM_PYTHON = Path("/Library/Frameworks/Python.framework/Versions/3.12/bin/python3")
VENV = PROJECT / ".venv-haddock3"
VENV_PYTHON = VENV / "bin" / "python"
PYMOL_APP = Path("/Applications/PyMOL.app")
PYMOL_BIN = PYMOL_APP / "Contents" / "MacOS" / "PyMOL"
LOCK_SYSTEM = PROJECT / "requirements-system-python-3.12.lock.txt"
LOCK_VENV = PROJECT / "requirements-haddock3.lock.txt"

HOST_ARCH = platform.machine()

# The load-bearing pins. Full lists live in the two .lock.txt files; these are the ones whose
# version changes a number or a launch path.
SYSTEM_PINS = {
    "rdkit": "2026.3.4", "scikit-learn": "1.9.0", "numpy": "2.5.1", "pandas": "2.3.3",
    "vina": "1.2.7", "meeko": "0.7.1", "requests": "2.34.2", "biopython": "1.87",
    "scipy": "1.18.0", "openpyxl": "3.1.5", "pdb-tools": "2.7.0",
}
VENV_PINS = {
    "haddock3": "2026.7.0", "pdb-tools": "2.7.0", "rdkit": "2026.3.4", "numpy": "2.5.1",
    "pandas": "2.3.3", "scipy": "1.18.0", "freesasa": "2.2.1", "gdock": "0.5.0",
    "prodigy-lig": "1.1.4", "prodigy-prot": "2.3.0",
}

# PROTOCOL_LOCK.md: 08 refuses to start without the first three; 09 needs the reference pair.
REQUIRED_INPUTS = [
    "ambig_FIXED.tbl",
    "crbn_active_FIXED.txt",
    "ppil4_active_FIXED.txt",
    "CRBN_receptor_thalidomide_Ryan.pdb",
    "PPIL4_chainB.pdb",
    "REFERENCE_FPFT.json",
]

# --------------------------------------------------------------------------- check bookkeeping

RESULTS: list[tuple[str, str, str]] = []


def check(name: str, status: str, detail: str = "") -> None:
    RESULTS.append((name, status, detail))


def run(cmd: list[str], timeout: int = 60) -> tuple[int, str, str]:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout.strip(), p.stderr.strip()
    except FileNotFoundError:
        return 127, "", "not found"
    except subprocess.TimeoutExpired:
        return 124, "", f"timed out after {timeout}s"
    except OSError as exc:
        return 1, "", str(exc)


def probe_pins(python: Path, pins: dict[str, str]) -> tuple[str, str]:
    """Return (status, detail) for a set of expected package versions in one interpreter."""
    code = (
        "import importlib.metadata as m, json\n"
        f"names = {sorted(pins)}\n"
        "out = {}\n"
        "for n in names:\n"
        "    try: out[n] = m.version(n)\n"
        "    except Exception: out[n] = None\n"
        "print(json.dumps(out))\n"
    )
    rc, out, err = run([str(python), "-c", code])
    if rc != 0:
        return "FAIL", f"cannot run {python}: {err.splitlines()[-1] if err else rc}"
    found = json.loads(out)
    missing = [n for n, v in pins.items() if found.get(n) is None]
    drifted = [f"{n} {found[n]} != {v}" for n, v in pins.items()
               if found.get(n) is not None and found[n] != v]
    if missing or drifted:
        bits = []
        if missing:
            bits.append("absent: " + ", ".join(missing))
        if drifted:
            bits.append("drifted: " + "; ".join(drifted))
        return "FAIL", " | ".join(bits)
    return "PASS", f"{len(pins)} pins match"


def stale_shebangs(bindir: Path) -> list[tuple[str, str]]:
    """Console scripts whose shebang does not point at this venv."""
    stale: list[tuple[str, str]] = []
    if not bindir.is_dir():
        return stale
    for f in sorted(bindir.iterdir()):
        if not f.is_file():
            continue
        try:
            first = f.open("rb").readline(300).decode("utf-8", "replace").strip()
        except OSError:
            continue
        if first.startswith("#!") and "python" in first:
            target = first[2:].strip()
            if str(VENV) not in target:
                stale.append((f.name, target))
    return stale


def lock_package_count(path: Path) -> int:
    if not path.exists():
        return -1
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines()
               if line.strip() and not line.startswith("#") and "==" in line)


# --------------------------------------------------------------------------- 1. host

check("host OS", "PASS", f"macOS {platform.mac_ver()[0]} ({platform.release()})")
check("host architecture", "PASS" if HOST_ARCH == "arm64" else "WARN",
      f"{HOST_ARCH}" + ("" if HOST_ARCH == "arm64" else " -- lock captured on arm64"))

rc, _, _ = run(["/usr/bin/arch", "-x86_64", "/usr/bin/true"])
check("Rosetta 2 (needed only for x86_64 PyMOL)", "PASS" if rc == 0 else "WARN",
      "available" if rc == 0 else "not installed -- x86_64 binaries cannot run")

# --------------------------------------------------------------------------- 2. interpreters

for label, exe in (("system Python", SYSTEM_PYTHON), ("venv Python", VENV_PYTHON)):
    if not exe.exists():
        check(f"{label} present", "FAIL", f"missing: {exe}")
        continue
    rc, out, err = run([str(exe), "-c", "import platform,sys; print(platform.python_version()); print(sys.executable)"])
    if rc != 0:
        check(f"{label} usable", "FAIL", err.splitlines()[-1] if err else str(rc))
        continue
    version, real = (out.splitlines() + ["", ""])[:2]
    check(f"{label} usable", "PASS" if version == "3.12.0" else "WARN",
          f"CPython {version}" + (f" at {real}" if label == "venv Python" else ""))

# --------------------------------------------------------------------------- 3. package pins

if SYSTEM_PYTHON.exists():
    status, detail = probe_pins(SYSTEM_PYTHON, SYSTEM_PINS)
    check("system Python pins", status, detail)
if VENV_PYTHON.exists():
    status, detail = probe_pins(VENV_PYTHON, VENV_PINS)
    check("venv pins", status, detail)

# --------------------------------------------------------------------------- 4. venv console scripts

stale = stale_shebangs(VENV / "bin")
if not (VENV / "bin").is_dir():
    check("venv console scripts", "WARN", f"no {VENV / 'bin'} -- venv not created yet?")
elif stale:
    old_venv = stale[0][1].rsplit("/bin/", 1)[0]
    check("venv console scripts", "FAIL",
          f"{len(stale)} point at a stale path, e.g. {stale[0][0]} -> {stale[0][1]}"
          f" | fix: sed -i '' 's|{old_venv}/|{VENV}/|' '{VENV}/bin/'*")
else:
    check("venv console scripts", "PASS", f"all point at {VENV}")

# --------------------------------------------------------------------------- 5. the CLI really runs

haddock3 = VENV / "bin" / "haddock3"
if haddock3.exists():
    rc, out, err = run([str(haddock3), "--version"])
    if rc == 0 and "2026.7.0" in out:
        check("haddock3 CLI", "PASS", out.replace("\n", " "))
    elif rc == 127:
        check("haddock3 CLI", "FAIL",
              "cannot exec: the interpreter named in its shebang does not exist "
              "(see the console-script check below)")
    else:
        check("haddock3 CLI", "FAIL",
              (err or out).splitlines()[-1] if (err or out) else f"exit {rc}")
else:
    check("haddock3 CLI", "FAIL", f"missing: {haddock3}")

# --------------------------------------------------------------------------- 6. PRODRG for this host

prodrg_dir = VENV / "lib" / "python3.12" / "site-packages" / "haddock" / "prodrg"
if prodrg_dir.is_dir():
    wanted = {"arm64": "prodrg_arm64-darwin", "x86_64": "prodrg_x86_64-darwin"}.get(HOST_ARCH)
    if wanted and (prodrg_dir / wanted).exists():
        check("bundled PRODRG", "PASS", f"{wanted} present")
    else:
        check("bundled PRODRG", "WARN", f"no build for {HOST_ARCH}; have: "
              + ", ".join(sorted(p.name for p in prodrg_dir.glob("prodrg_*"))))
else:
    check("bundled PRODRG", "WARN", "haddock3 not installed in the venv")

# --------------------------------------------------------------------------- 7. PyMOL

if not PYMOL_BIN.exists():
    check("PyMOL.app", "FAIL", f"missing: {PYMOL_BIN}")
else:
    try:
        info = plistlib.load((PYMOL_APP / "Contents" / "Info.plist").open("rb"))
        version = info.get("CFBundleShortVersionString", "?")
    except Exception as exc:                                    # noqa: BLE001
        version = f"? ({exc})"
    rc, archs, _ = run(["lipo", "-archs", str(PYMOL_BIN)])
    archs = archs or "?"
    rc2, out2, err2 = run([str(PYMOL_BIN), "-cq", "print('pymol-ok')"], timeout=120)
    if rc2 == 0:
        check("PyMOL.app launchable", "PASS", f"version {version}, arch {archs}")
    else:
        check("PyMOL.app launchable", "FAIL",
              f"version {version}, arch {archs} on {HOST_ARCH} host -- "
              f"{(err2 or out2).splitlines()[-1] if (err2 or out2) else f'exit {rc2}'}"
              " | stages 07/13/15/16 cannot run; needs an arm64 PyMOL.app")

# --------------------------------------------------------------------------- 8. protocol inputs

missing = [n for n in REQUIRED_INPUTS if not (PROJECT / n).exists()]
check("protocol input files", "FAIL" if missing else "PASS",
      "missing: " + ", ".join(missing) if missing else f"{len(REQUIRED_INPUTS)} present")

# --------------------------------------------------------------------------- 9. lock coverage

for label, lock, exe in (("system", LOCK_SYSTEM, SYSTEM_PYTHON), ("venv", LOCK_VENV, VENV_PYTHON)):
    if not lock.exists():
        check(f"{label} lock file", "FAIL", f"missing: {lock.name}")
        continue
    rc, out, err = run([str(exe), "-m", "pip", "list", "--format=freeze"])
    if rc != 0:
        check(f"{label} lock file", "WARN", f"cannot count installed packages: {err}")
        continue
    installed = len([ln for ln in out.splitlines() if "==" in ln])
    pinned = lock_package_count(lock)
    check(f"{label} lock covers installed packages", "PASS" if pinned == installed else "WARN",
          f"{pinned} pinned, {installed} installed")

# --------------------------------------------------------------------------- report

width = max(len(n) for n, _, _ in RESULTS)
order = {"FAIL": 0, "WARN": 1, "PASS": 2}
print()
for name, status, detail in sorted(RESULTS, key=lambda r: (order[r[1]], r[0])):
    print(f"  {status:<4}  {name:<{width}}  {detail}")
counts = {s: sum(1 for _, st, _ in RESULTS if st == s) for s in ("PASS", "WARN", "FAIL")}
print(f"\n  {counts['PASS']} passed, {counts['WARN']} warning(s), {counts['FAIL']} failed"
      f"   -- ENVIRONMENT_LOCK.md (captured 2026-10-11)")
if counts["FAIL"]:
    print("  FIX THE FAILURES BEFORE DOCKING: a broken interpreter fails mid-campaign, not at launch.\n")
    sys.exit(1)
print()
sys.exit(0)
