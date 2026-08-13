#!/usr/bin/env python3
"""Run the locked-environment audit with narrow, documented Django 4.2 triage."""

from __future__ import annotations

import subprocess
import sys

# Django 4.2.30 is the final required 4.2-series baseline for issue #15. These
# advisories have no 4.2 patch. They are low/medium, and their affected APIs are
# absent from this repository; see docs/operations/DEPENDENCY_AUDIT_TRIAGE.md.
TRIAGED_DJANGO_42_ADVISORIES = (
    "GHSA-923m-gv2p-w5qp",
    "GHSA-h7pc-vwp9-298g",
    "GHSA-8cjm-8mp7-r2xf",
    "GHSA-3h9f-r86x-qvjx",
    "GHSA-crhf-3pfg-w68w",
    "GHSA-8qcx-xf44-272x",
)

# pytest-playwright requires pytest <9. The advisory requires a malicious local
# user on a shared UNIX host; tests run in isolated developer/CI environments.
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
    for advisory in TRIAGED_DJANGO_42_ADVISORIES + TRIAGED_TEST_TOOL_ADVISORIES:
        command.extend(("--ignore-vuln", advisory))
    return subprocess.run(command, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
