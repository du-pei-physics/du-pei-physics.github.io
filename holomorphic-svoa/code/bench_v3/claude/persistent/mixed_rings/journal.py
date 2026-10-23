"""Crash consistency for the benchmark-v3 workers (Codex review of 19 Oct 2026, item 6), 19 Oct 2026.

The state JSON, rowspaces.txt and hit_index.jsonl are separate writes, so an abrupt interruption can leave the counters
in state.json behind the journals.  On every start the workers call reconcile(): the hit counters are recomputed from
hit_index.jsonl (the durable journal: one line per exact hit, with its verification outcome and rowspace hash), and the
distinct count from rowspaces.txt; any disagreement with the checkpointed counters is returned so that the worker
records it in runs.jsonl.  A missing "verified" field counts as not verified; a malformed line is counted and reported.
"""
import json
from pathlib import Path


def reconcile(folder, st):
    folder = Path(folder)
    hi = folder / "hit_index.jsonl"
    rs = folder / "rowspaces.txt"
    j = {"hits": {"nega": 0, "cyc": 0}, "verified": 0, "verify_fail": 0, "cand_hit": 0, "malformed": 0}
    cand, hashes_in_index = set(), set()
    if hi.exists():
        for line in hi.read_text().splitlines():
            try:
                r = json.loads(line)
                kind = r["kind"]
                j["hits"][kind] += 1
            except Exception:
                j["malformed"] += 1
                continue
            if r.get("verified") is True:
                j["verified"] += 1
                if r.get("rowspace_sha256"):
                    hashes_in_index.add(r["rowspace_sha256"])
            else:
                j["verify_fail"] += 1
            cand.add((r.get("run"), r.get("epoch"), tuple(r.get("tuple") or ())))
    j["cand_hit"] = len(cand)
    rows = set(x for x in rs.read_text().split() if x) if rs.exists() else set()
    diff = {}
    for key in ("verified", "verify_fail", "cand_hit"):
        if st.get(key) != j[key]:
            diff[key] = {"state": st.get(key), "journal": j[key]}
        st[key] = j[key]
    if st.get("hits") != j["hits"]:
        diff["hits"] = {"state": st.get("hits"), "journal": j["hits"]}
    st["hits"] = j["hits"]
    if st.get("distinct") != len(rows):
        diff["distinct"] = {"state": st.get("distinct"), "rowspaces": len(rows)}
    st["distinct"] = len(rows)
    if rows - hashes_in_index:
        diff["rowspaces_without_index_line"] = len(rows - hashes_in_index)
    if j["malformed"]:
        diff["malformed_index_lines"] = j["malformed"]
    return diff, rows
