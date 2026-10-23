# Copied from claude/persistent/common.py of the scft-code-search workspace for the holomorphic-svoa note;
# only import/output paths were changed (see tools/build_code_folder.py).
"""Shared infrastructure for Claude's persistent workers (see README.md)."""
import json, os, shutil, time
from pathlib import Path

MIN_FREE_GIB = 10


class Worker:
    def __init__(self, folder, name):
        self.folder = Path(folder); self.folder.mkdir(parents=True, exist_ok=True)
        self.name = name
        self.state_path = self.folder / "state.json"
        self.stop_path = self.folder / "stop.flag"
        self.state = {"name": name, "pid": os.getpid(), "started": time.strftime("%Y-%m-%d %H:%M:%S"),
                      "attempts": 0, "hits": 0, "last_update": None, "status": "running"}
        if self.state_path.exists():
            try:
                old = json.loads(self.state_path.read_text())
                self.state["attempts"] = old.get("attempts", 0); self.state["hits"] = old.get("hits", 0)
                self.state["previous_runs"] = old.get("previous_runs", 0) + 1
            except Exception:
                pass
        self._last_write = 0
        self.save(force=True)

    def save(self, force=False, **extra):
        self.state.update(extra)
        if force or time.time() - self._last_write > 20:
            self.state["last_update"] = time.strftime("%Y-%m-%d %H:%M:%S")
            tmp = self.state_path.with_suffix(".tmp")
            tmp.write_text(json.dumps(self.state, indent=1)); os.replace(tmp, self.state_path)
            self._last_write = time.time()

    def should_stop(self):
        return self.stop_path.exists()

    def wait_for_disk(self):
        while shutil.disk_usage(self.folder).free < MIN_FREE_GIB * 2**30:
            self.save(force=True, status="paused: low disk")
            if self.should_stop():
                return
            time.sleep(300)
        if self.state.get("status", "").startswith("paused"):
            self.save(force=True, status="running")

    def append(self, fname, record):
        with open(self.folder / fname, "a") as f:
            f.write(json.dumps(record) + "\n")

    def finish(self, status="finished"):
        self.save(force=True, status=status)
