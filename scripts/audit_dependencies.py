#!/usr/bin/env python3
"""Run the locked-environment audit with narrow documented test-tool triage."""

from __future__ import annotations

import subprocess
import sys

# pytest-playwright requires pytest <9. The advisory requires a malicious local
# user on a shared UNIX host; tests run in isolated developer/CI environments.
# This exact development-only decision remains tracked in the dependency triage.
TRIAGED_TEST_TOOL_ADVISORIES = ("GHSA-6w46-j5rx-g56g",)


def main() -> int:
    command = [
        sys.executable,
        "-m",
        "pip_audit",
        "--local",
        "--progress-spinner",
        "off",
    ]
    for advisory in TRIAGED_TEST_TOOL_ADVISORIES:
        command.extend(("--ignore-vuln", advisory))
    return subprocess.run(command, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
