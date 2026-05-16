from __future__ import annotations

import re
from pathlib import Path
from typing import Any


SENSITIVE_PATH_PARTS = {
    "local_private",
    "media",
}
SENSITIVE_SUFFIXES = {
    ".sqlite3",
    ".db",
    ".pdf",
    ".png",
    ".jpg",
    ".jpeg",
    ".tif",
    ".tiff",
}
ALLOWED_PSEUDONYM_PATTERNS = (
    re.compile(r"Research Patient\d*", re.IGNORECASE),
    re.compile(r"MM_RESEARCH_[A-Z0-9_]+", re.IGNORECASE),
)
_DIRECT_IDENTIFIER_TERMS = (
    "Ros" + "sana",
    "Agu" + "eci",
)
PATIENT_NAME_PATTERNS = tuple(
    re.compile(rf"\b{re.escape(term)}\b", re.IGNORECASE) for term in _DIRECT_IDENTIFIER_TERMS
)
TAX_CODE_PATTERN = re.compile(r"\b[A-Z]{6}\d{2}[A-Z]\d{2}[A-Z]\d{3}[A-Z]\b")
PHONE_PATTERN = re.compile(r"(?<!\w)(?:\+?39[\s.-]?)?(?:\d[\s.-]?){8,12}\d(?!\w)")
SENSITIVE_TEXT_MARKERS = (
    "Codice" + " Fiscale",
    "Data" + " di nascita",
    "birth_date",
    "medical_record_number",
)


def scan_text_for_sensitive_markers(text: str) -> list[dict[str, str]]:
    findings: list[dict[str, str]] = []
    scrubbed = text
    for allowed in ALLOWED_PSEUDONYM_PATTERNS:
        scrubbed = allowed.sub("", scrubbed)
    for pattern in PATIENT_NAME_PATTERNS:
        if pattern.search(scrubbed):
            findings.append({"status": "fail", "check": "patient_name", "detail": pattern.pattern})
    if TAX_CODE_PATTERN.search(scrubbed):
        findings.append({"status": "fail", "check": "tax_code_pattern", "detail": "tax-code-like string"})
    if PHONE_PATTERN.search(scrubbed):
        findings.append({"status": "fail", "check": "phone_pattern", "detail": "phone-like string"})
    for marker in SENSITIVE_TEXT_MARKERS:
        if marker.lower() in scrubbed.lower():
            findings.append({"status": "fail", "check": "sensitive_marker", "detail": marker})
    return findings


def is_sensitive_path(path: str | Path) -> bool:
    value = Path(str(path))
    parts = {part.lower() for part in value.parts}
    if parts & SENSITIVE_PATH_PARTS:
        return True
    if value.name in {"db.sqlite3", "db.sqlite3-journal", ".env"}:
        return True
    return value.suffix.lower() in SENSITIVE_SUFFIXES


def check_path_payload(path: str | Path, text: str | None = None) -> list[dict[str, Any]]:
    findings: list[dict[str, Any]] = []
    if is_sensitive_path(path):
        findings.append({"status": "fail", "check": "sensitive_path", "object_id": str(path), "next_action": "Do not stage or push this path."})
    if text:
        findings.extend(scan_text_for_sensitive_markers(text))
    return findings