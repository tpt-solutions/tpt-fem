#!/usr/bin/env python3
"""Check that the crates excluded from the Cargo workspace (py/capi/wasm)
pin the same `tpt-fem*` versions as the root `[workspace.dependencies]`.

Cargo does not validate these (the crates are not workspace members), so a
version bump in the workspace can silently leave them building against a stale
pin. Exits non-zero and lists every mismatch.
"""
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXCLUDED = ["tpt-fem-py", "tpt-fem-capi", "tpt-fem-wasm"]


def main() -> int:
    ws = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))
    pins = {
        name: spec["version"]
        for name, spec in ws["workspace"]["dependencies"].items()
        if name.startswith("tpt-fem") and isinstance(spec, dict) and "version" in spec
    }
    problems = []
    for crate in EXCLUDED:
        manifest = tomllib.loads(
            (ROOT / "crates" / crate / "Cargo.toml").read_text(encoding="utf-8")
        )
        for name, spec in manifest.get("dependencies", {}).items():
            if not name.startswith("tpt-fem"):
                continue
            version = spec.get("version") if isinstance(spec, dict) else spec
            if name not in pins:
                problems.append(f"{crate}: depends on unknown workspace crate {name}")
            elif version != pins[name]:
                problems.append(
                    f"{crate}: {name} pinned to {version!r}, workspace has {pins[name]!r}"
                )
    for p in problems:
        print("MISMATCH:", p)
    if not problems:
        print("excluded-crate pins match the workspace")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
