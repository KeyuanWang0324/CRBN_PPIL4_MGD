# ENVIRONMENT_LOCK — the three interpreters, pinned

Captured **2026-10-11** on this machine. This is the environment counterpart to
[`PROTOCOL_LOCK.md`](PROTOCOL_LOCK.md): that file pins the *docking protocol*, this one pins the
*software that runs it*. `PROTOCOL_LOCK.md` §3 pins exactly one version ("HADDOCK3 2026.7.0");
until now nothing recorded rdkit, scikit-learn, Vina, Meeko, or PyMOL — the four libraries every
number in `08_*` inherits.

**Rule 0 of the protocol — "share the artifact, not a description" — applies here too.** The two
`.lock.txt` files are the artifact; this document is the description.

## 0. Status on this machine right now

| Component | State | Evidence |
|---|---|---|
| System Python 3.12 (stages 00–04, 11, 12, 14) | **OK** | 81 packages import; `rdkit 2026.3.4`, `vina 1.2.7` |
| `.venv-haddock3` packages (05, 06, 08, 09) | **OK** — HADDOCK3 2026.7.0, all bundled binaries are **arm64** | `prodrg_arm64-darwin`, `contact_fcc`, `fast-rmsdmatrix`, `haddock-restraints` → `Mach-O arm64` |
| `.venv-haddock3` **console scripts** | **BROKEN** — 71 of the 78 files in `bin/` are console scripts whose shebang points at the pre-move path | `head -1 .venv-haddock3/bin/haddock3` → `#!/Users/ryanwang/Desktop/ISEF/.venv-haddock3/bin/python` |
| PyMOL.app (stages 07, 13, 15, 16) | **BROKEN** — x86_64 app on an arm64 host, and Rosetta is absent | `PyMOL -cq` → `Bad CPU type in executable`; `arch -x86_64 true` → same |
| `pymol-open-source` in the system Python | **BROKEN** (as README documents) | `import pymol` → `Library not loaded: @rpath/libpng16.16.dylib` from `/Users/Martin/.local/…/pymol_312/lib` |

Three of the four breakages are repairable in minutes (§6). The PyMOL one blocks re-running the
structure-figure and viewer stages — note that does **not** endanger the committed deliverables:
`ternary_structures/*.pse`, `CRBN_Project_site/glb/*` and `img/*` already exist; only
*regeneration* is blocked.

## 1. Host

| | |
|---|---|
| OS | macOS **27.0.1** (build 26A434) |
| Architecture | **arm64** |
| Rosetta 2 | **not available** — `arch -x86_64 /usr/bin/true` fails with `Bad CPU type in executable` |
| Sync caveat | the project lives under `~/Desktop/校外/`, i.e. **inside iCloud sync** — see `PROTOCOL_LOCK.md` §3a. Any `run_dir` must stay outside it. |

## 2. The three interpreters

PyMOL is not a library here, it is its own interpreter — that is why the pipeline needs three.

| Stage | Interpreter | Absolute path | Version | Lock file |
|---|---|---|---|---|
| 00–04, 11, 12, 14 | system Python | `/Library/Frameworks/Python.framework/Versions/3.12/bin/python3` | CPython **3.12.0** (pip 23.2.1) | [`requirements-system-python-3.12.lock.txt`](requirements-system-python-3.12.lock.txt) — 81 packages |
| 05, 06, 08, 09 | HADDOCK3 venv | `.venv-haddock3/bin/python` → the same 3.12.0 framework build | CPython **3.12.0** (pip 26.1.2) | [`requirements-haddock3.lock.txt`](requirements-haddock3.lock.txt) — 28 packages |
| 07, 13, 15, 16 | PyMOL's bundled Python | `/Applications/PyMOL.app/Contents/MacOS/PyMOL` | PyMOL **3.1.6.1** (bundle is x86_64) | not pip-managed; version from `Info.plist` |

Most scripts detect the wrong interpreter and re-exec; `08` additionally prepends
`.venv-haddock3/bin` to `PATH` and calls the `haddock3` console script (`08:93-98`, `08:632`) —
which is exactly what §6a breaks.

## 3. Package pins that matter

Only the load-bearing pins are listed; the `.lock.txt` files carry the full `pip freeze`.

**System Python** (81 packages) — the numbers in `02/03` and `04` depend on these:

| Package | Version | Used by |
|---|---|---|
| `rdkit` | 2026.3.4 | 01, 11, 12, 14 (library, descriptors, canonical SMILES) |
| `scikit-learn` | 1.9.0 | 02, 03 (RandomForestClassifier) |
| `numpy` / `pandas` | 2.5.1 / 2.3.3 | throughout |
| `vina` | **1.2.7** | 04 (Python bindings; there is **no** `vina` CLI on `PATH`) |
| `meeko` | 0.7.1 | 04 (PDBQT preparation) |
| `scipy` | 1.18.0 | statistics in the analysis scripts |
| `requests` | 2.34.2 | 02 (ChEMBL REST) |
| `biopython` | 1.87 | structure handling |
| `openpyxl` | 3.1.5 | xlsx workbooks |
| `pdb-tools` | 2.7.0 | chain surgery |
| `pymol-open-source` | 3.2.0a0 | **installed but not importable** — do not rely on it |

**HADDOCK3 venv** (28 packages):

| Package | Version | Note |
|---|---|---|
| `haddock3` | **2026.7.0** | the one version `PROTOCOL_LOCK.md` §3 pins — do not upgrade mid-campaign |
| `pdb-tools` | 2.7.0 | console scripts affected by §6a |
| `freesasa` | 2.2.1 | SASA/BSA in the structure axis |
| `gdock` | 0.5.0 | HADDOCK dependency |
| `prodigy-lig` / `prodigy-prot` | 1.1.4 / 2.3.0 | contact-based scoring |
| `rdkit` | 2026.3.4 | matches the system Python |
| `numpy` / `pandas` / `scipy` | 2.5.1 / 2.3.3 / 1.18.0 | matches |
| `jsonpickle` / `toml` / `PyYAML` / `psutil` | 4.1.2 / 0.10.2 / 6.0.3 / 7.2.2 | HADDOCK3 internals |

> The two environments are deliberately **aligned** on `rdkit`, `numpy`, `pandas`, `scipy`,
> `biopython` and `pdb-tools`. Keep them aligned: `11` re-derives columns that `08` produced, and
> a cross-environment descriptor difference would be silently attributed to the science.

## 4. Bundled and external tools

| Tool | Where | Architecture check |
|---|---|---|
| HADDOCK3 CLI | `.venv-haddock3/bin/haddock3` | python script (shebang — see §6a) |
| PRODRG | `…/site-packages/haddock/prodrg/prodrg_arm64-darwin` | **arm64 — correct for this host** |
| `contact_fcc`, `fast-rmsdmatrix`, `haddock-restraints` | `…/haddock/bin/` | arm64 / universal ✓ |
| AutoDock Vina | Python module only (`vina==1.2.7`); no CLI | n/a |
| PyMOL | `/Applications/PyMOL.app` (3.1.6.1) | **x86_64 — wrong for this host** |

HADDOCK3 ships four PRODRG builds (`prodrg_arm64-darwin`, `prodrg_aarch64-linux`,
`prodrg_x86_64-linux`, plus `prodrg.param`) and selects by platform; on this machine it picks the
arm64 one, so **the docking stack is native and unaffected by the Rosetta problem.**

## 5. Recreating an environment from scratch

```bash
# --- system Python (stages 00–04, 11, 12, 14) ---
/Library/Frameworks/Python.framework/Versions/3.12/bin/python3 -m pip install \
    -r requirements-system-python-3.12.lock.txt

# --- HADDOCK3 venv (stages 05, 06, 08, 09) ---
/Library/Frameworks/Python.framework/Versions/3.12/bin/python3.12 -m venv .venv-haddock3
.venv-haddock3/bin/python -m pip install -r requirements-haddock3.lock.txt
.venv-haddock3/bin/haddock3 --version      # must print 2026.7.0

# --- PyMOL (stages 07, 13, 15, 16) ---
# Not pip-installable on this machine. Needs an arm64 PyMOL.app — see §6b.
```

After either install, re-capture the pin:

```bash
/Library/Frameworks/Python.framework/Versions/3.12/bin/python3 -m pip freeze > requirements-system-python-3.12.lock.txt
.venv-haddock3/bin/python -m pip freeze > requirements-haddock3.lock.txt
```

## 6. Known breakages, with the fix verified on this machine

### 6a. `.venv-haddock3` console scripts — 71 of 78 have a stale shebang

The project directory moved from `~/Desktop/ISEF` to `~/Desktop/校外/ISEF`. A venv hard-codes its
own path at creation time, in the shebang of every console script, so `haddock3` and all of
`pdb-tools` now fail with:

```
bad interpreter: No such file or directory
```

74 files under `bin/` mention the old path; 71 of them are python console scripts and the other
three are the `activate` / `activate.csh` / `activate.fish` shell scripts, which carry it as a
comment. `pyvenv.cfg` still records the old path too
(`command = … -m venv /Users/ryanwang/Desktop/ISEF/.venv-haddock3`).

**Fix (seconds, verified).** Rewrite the shebangs; nothing else inside the venv needs touching —
a scan found **0** stale references outside `bin/` excluding `.pyc` bytecode caches:

```bash
sed -i '' 's|/Users/ryanwang/Desktop/ISEF/|/Users/ryanwang/Desktop/校外/ISEF/|' .venv-haddock3/bin/*
.venv-haddock3/bin/haddock3 --version     # -> haddock3 - 2026.7.0
```

`00_verify_environment_lock_(Ryan).py` prints this exact command, with the old and new paths
filled in from whatever it finds. Verified before recommending: a repaired copy run from `/tmp`
printed `haddock3 - 2026.7.0` (rc 0). The alternative is §5's "recreate the venv", which is
slower and can pull different versions — **prefer the sed.** `os.path.join(SCRIPT_DIR,
".venv-haddock3", "bin")` means the scripts will find the venv automatically once its scripts
work again.

### 6b. PyMOL.app is x86_64 on an arm64 host, and Rosetta is gone

```
$ file /Applications/PyMOL.app/Contents/MacOS/PyMOL
Mach-O 64-bit executable x86_64
$ /Applications/PyMOL.app/Contents/MacOS/PyMOL -cq "print(1)"
Bad CPU type in executable
$ /usr/bin/arch -x86_64 /usr/bin/true
Bad CPU type in executable          # Rosetta 2 is not installed/available
```

`MacPyMOL`, `PyMOL` and the bundled `python` are all x86_64; there is no arm64 slice. So
**07, 13, 15 and 16 cannot run on this machine**, and the `.pse` files, the site's `glb/` models
and the static renders cannot be regenerated until PyMOL is replaced with an **arm64 build**
(or Rosetta is restored). Every committed output of those stages is intact — this is a
regeneration problem, not a data-loss problem. The README's claim that "the installed PyMOL.app
works headlessly" was true when written and is **no longer true**; it should be re-checked
before anyone reruns 07/13/15/16.

### 6c. `pymol-open-source` in the system Python

`pymol-open-source==3.2.0a0` is installed but unimportable: the compiled `_cmd` extension
hard-codes an rpath into the wheel-builder's mamba env
(`/Users/Martin/.local/share/mamba/envs/pymol_312/lib/libpng16.16.dylib`). Reinstalling does not
fix it, and `pip install` a second time will not either — this is what README already warns
about. Keep it out of the lock-driven workflow; it is recorded in the lock file only so a
`pip install -r` reproduces the machine exactly.

## 7. Verification

```bash
/Library/Frameworks/Python.framework/Versions/3.12/bin/python3 "00_verify_environment_lock_(Ryan).py"
```

[`00_verify_environment_lock_(Ryan).py`](<00_verify_environment_lock_(Ryan).py>) checks, in order:
host arch and OS; each interpreter's path and exact version; the load-bearing package pins in
both environments; the venv shebang integrity (§6a) and whether `haddock3 --version` really runs;
the PyMOL bundle's architecture and launchability (§6b); the presence of the arm64 PRODRG build;
the protocol's required input files (`ambig_FIXED.tbl`, `*_active_FIXED.txt`,
`CRBN_receptor_thalidomide_Ryan.pdb`, `PPIL4_chainB.pdb`, the 9DWV reference, `REFERENCE_FPFT.json`);
and that the two lock files cover the installed package sets.

Exit **0** = every check passed. Exit **1** = at least one FAIL (warnings alone do not fail).
Run it before a docking campaign, and after any folder move, Python upgrade or `pip install` —
it is the check that would have caught §6a on the day the folder moved, instead of weeks later.

## 8. Honest limits of this lock

- **It pins this machine, not a portable spec.** The system Python is a python.org framework
  build with a global `pip` (23.2.1); there is no `pyproject.toml`, no hash-pinned lock
  (`pip-compile`/`uv`), and no CI. Reproducing this on Linux would need a new capture.
- **No hashes.** `==` pins are exact-version, not byte-exact artifacts. If a rebuild ever
  matters, add `pip-compile --generate-hashes`.
- **PyMOL is not captured by pip at all** — only by its `Info.plist` version, which cannot be
  installed from this file.
- **The `dockq_vs_9dwv` blind spot is unrelated to this lock**: no second raw run is committed,
  so "byte-identical" reproduction cannot be verified from the repository even with a perfect
  environment (§`PROTOCOL_LOCK.md` §5 specifies tolerance, not byte equality).
- macOS itself is not pinned: **27.0.1 / arm64** is recorded, and the arm64 requirement is
  load-bearing given §6b.
