#!/usr/bin/env python3
"""Run the locked-environment audit with narrow, expiring advisory triage."""

from __future__ import annotations

import subprocess
import sys
from datetime import date

# Django 4.2 reached end of extended support on 2026-04-07. These exact
# exceptions are temporary risk decisions only; they do not make the framework
# supported or authorize shared/production promotion. See issues #69 and #70
# and docs/operations/DEPENDENCY_AUDIT_TRIAGE.md.
TEMPORARY_DJANGO_42_ADVISORIES = (
    "GHSA-923m-gv2p-w5qp",
    "GHSA-h7pc-vwp9-298g",
    "GHSA-8cjm-8mp7-r2xf",
    "GHSA-3h9f-r86x-qvjx",
    "GHSA-crhf-3pfg-w68w",
    "GHSA-8qcx-xf44-272x",
    "PYSEC-2026-3717",
)
TEMPORARY_DJANGO_42_EXCEPTION_EXPIRES = date(2026, 9, 30)

# pytest-playwright requires pytest <9. The advisory requires a malicious local
# user on a shared UNIX host; tests run in isolated developer/CI environments.
TRIAGED_TEST_TOOL_ADVISORIES = ("GHSA-6w46-j5rx-g56g",)


def _temporary_django_exception_is_current(today: date) -> bool:
    return today <= TEMPORARY_DJANGO_42_EXCEPTION_EXPIRES


def main() -> int:
    today = date.today()
    if not _temporary_django_exception_is_current(today):
        print(
            "ERROR: temporary Django 4.2 advisory triage expired on "
            f"{TEMPORARY_DJANGO_42_EXCEPTION_EXPIRES.isoformat()}; complete issue #70 "
            "and remove the unsupported-framework exceptions.",
            file=sys.stderr,
        )
        return 2

    command = [
        sys.executable,
        "-m",
        "pip_audit",
        "--local",
        "--progress-spinner",
        "off",
    ]
    for advisory in TEMPORARY_DJANGO_42_ADVISORIES + TRIAGED_TEST_TOOL_ADVISORIES:
        command.extend(("--ignore-vuln", advisory))
    return subprocess.run(command, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
