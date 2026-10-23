"""Restore the worker sources of one recorded benchmark-v3 run into a scratch directory and print its replay command.

Usage (in code/bench_v3/):  python3 restore_run.py LANE RUN SCRATCH      e.g.  python3 restore_run.py bench3_k21b 3 /tmp/replay

LANE is a lane name (its run records are runs/LANE.jsonl), RUN a run number in that file, SCRATCH a new or empty directory.
For every source file whose SHA-256 the run recorded, the archived file with exactly that hash (SOURCES.json says which
one) is copied to its original workspace-relative path under SCRATCH, and the hash is checked again after copying.  A file
with a documented publication edit (SOURCES.json: publication_edit, packaged_sha256; a docstring phrase only) is checked
against its packaged hash instead, and reported.  The
dependencies whose version the runs did not record (SOURCES.json, unrecorded_dependencies) are copied at their archived
version, with a warning.  SCRATCH/RESTORED.json lists what was restored.  The script then prints the replay command: the
recorded argv, with <workspace> mapped to SCRATCH, to be started in SCRATCH/claude/persistent/mixed_rings; the worker
writes its lane folder next to mixed_rings/ inside SCRATCH.

Scope.  A run draws its pools and candidates from numpy's default_rng with the run's recorded seed, so a replay from an
empty lane folder repeats the run's sequence of candidates; the epoch numbering restarts at 1 for a run that resumed a lane,
the --hits-max budget of stored full records starts afresh, and where the replay stops depends on --max-minutes and the
machine.  The RNG state before each pool is in the lane's epochs.jsonl (audit bundle in ../../data/), for comparison.
Not reproduced: the run-time versions of the unrecorded dependencies, the Python and numpy builds (the run record gives
their versions; install them yourself, with sympy) and the rest of the historical environment.
"""
import hashlib, json, shlex, shutil, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(argv):
    if len(argv) != 4 or not argv[2].lstrip("-").isdigit():
        raise SystemExit(__doc__.split("\n\n")[1])
    lane, run_no, scratch = argv[1], int(argv[2]), Path(argv[3]).resolve()
    sources = json.loads((HERE / "SOURCES.json").read_text())
    runs_file = HERE / "runs" / f"{lane}.jsonl"
    if not runs_file.exists():
        raise SystemExit(f"no run records for lane {lane!r} (runs/ has: "
                         + ", ".join(sorted(p.stem for p in (HERE / "runs").glob("*.jsonl"))) + ")")
    recs = [r for r in map(json.loads, filter(str.strip, runs_file.read_text().splitlines())) if r.get("run") == run_no]
    if len(recs) != 1:
        raise SystemExit(f"{lane}: {len(recs)} run records with run number {run_no}")
    rec = recs[0]
    if scratch.exists() and any(scratch.iterdir()):
        raise SystemExit(f"{scratch} is not empty; give a new or empty directory")
    restored, unrecorded, edited = [], [], []
    for name, h in sorted((rec.get("sources") or {}).items()):
        entry = sources["recorded"].get(name, {}).get(h)
        if not entry or not entry.get("packaged"):
            raise SystemExit(f"{name} with SHA-256 {h[:12]}...: this version is not in the archive (see SOURCES.json)")
        src = HERE / entry["packaged"]
        expect = entry.get("packaged_sha256", h)             # differs from h only for a documented publication edit
        if sha(src) != expect:
            raise SystemExit(f"{entry['packaged']}: SHA-256 differs from the run record and SOURCES.json")
        dest = scratch / entry["workspace_path"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dest)
        if sha(dest) != expect:
            raise SystemExit(f"{dest}: SHA-256 differs after copying")
        restored.append({"path": entry["workspace_path"], "sha256": h, "from": entry["packaged"]})
        if expect != h:
            restored[-1].update({"packaged_sha256": expect, "publication_edit": entry.get("publication_edit")})
            edited.append(restored[-1])
    for rel, entry in sorted((sources.get("unrecorded_dependencies") or {}).items()):
        dest = scratch / rel
        if dest.exists():
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(HERE / entry["packaged"], dest)
        unrecorded.append({"path": rel, "sha256": sha(dest), "note": "archived version; its version at run time was not recorded"})
    cmd = [a.replace("<workspace>", str(scratch)) for a in rec.get("argv", [])]
    workdir = scratch / "claude/persistent/mixed_rings"
    (scratch / "RESTORED.json").write_text(json.dumps({
        "lane": lane, "run": run_no, "run_record": rec, "restored_recorded_sources": restored,
        "unrecorded_dependencies": unrecorded, "publication_edits": edited, "workdir": str(workdir),
        "command": ["python3"] + cmd}, indent=1) + "\n")
    print(f"{lane} run {run_no} ({rec.get('time')}): {len(restored)} recorded source files restored, each checked against its SHA-256:")
    for r in restored:
        print(f"  {r['path']}  {r['sha256'][:12]}")
    if edited:
        print(f"NOTE: {len(edited)} restored file(s) differ from the recorded SHA-256 by a documented publication edit:")
        for r in edited:
            print(f"  {r['path']}: {r['publication_edit']}")
    if unrecorded:
        print(f"WARNING: {len(unrecorded)} dependencies were copied at their archived version; the version they had when "
              "the run started was not recorded:")
        for r in unrecorded:
            print(f"  {r['path']}")
    try:
        import numpy
        here_np = numpy.__version__
    except ImportError:
        here_np = "not installed"
    print(f"recorded environment: Python {rec.get('python')}, numpy {rec.get('numpy')}; "
          f"this interpreter: Python {sys.version.split()[0]}, numpy {here_np}")
    print("replay (the lane folder is created in SCRATCH/claude/persistent/; adjust --max-minutes as needed):")
    print(f"  cd {shlex.quote(str(workdir))} && python3 " + " ".join(shlex.quote(a) for a in cmd))


if __name__ == "__main__":
    main(sys.argv)
