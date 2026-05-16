#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

python_bin="${PYTHON_BIN:-./venv/bin/python}"
if [[ ! -x "$python_bin" ]]; then
  echo "ERROR: Python interpreter not found or not executable: $python_bin" >&2
  exit 1
fi

echo "== Research safety check =="

required_ignored=("local_private" "media" "db.sqlite3" ".env")
for path in "${required_ignored[@]}"; do
  if ! git check-ignore -q "$path"; then
    echo "ERROR: expected ignored path is not ignored: $path" >&2
    exit 1
  fi
done

deleted_sensitive="$( { git diff --name-only --diff-filter=D; git diff --cached --name-only --diff-filter=D; } | sort -u || true)"
tracked_sensitive="$(git ls-files | grep -E '(^|/)(local_private|media|logs)(/|$)|(^|/)db\.sqlite3$|\.(sqlite3|db|log|pem|key|p12)$' || true)"
if [[ -n "$tracked_sensitive" && -n "$deleted_sensitive" ]]; then
  tracked_sensitive="$(comm -23 <(printf '%s\n' "$tracked_sensitive" | sort -u) <(printf '%s\n' "$deleted_sensitive" | sort -u) || true)"
fi
if [[ -n "$tracked_sensitive" ]]; then
  echo "ERROR: tracked sensitive paths detected:" >&2
  echo "$tracked_sensitive" >&2
  exit 1
fi

staged_files="$(git diff --cached --name-only --diff-filter=ACMRT || true)"
if [[ -n "$staged_files" ]]; then
  staged_sensitive="$(printf '%s\n' "$staged_files" | grep -E '(^|/)(local_private|media|logs)(/|$)|(^|/)db\.sqlite3$|\.(sqlite3|db|pem|key|p12|pdf|png|jpg|jpeg|tif|tiff)$' || true)"
  if [[ -n "$staged_sensitive" ]]; then
    echo "ERROR: staged sensitive paths detected:" >&2
    echo "$staged_sensitive" >&2
    exit 1
  fi

  staged_phi=""
  direct_name_one="Ros""sana"
  direct_name_two="Agu""eci"
  tax_marker="Codice"" Fiscale"
  birth_marker="Data"" di nascita"
  phi_regex="\\b(${direct_name_one}|${direct_name_two})\\b|${tax_marker}|${birth_marker}|\\b[A-Z]{6}[0-9]{2}[A-Z][0-9]{2}[A-Z][0-9]{3}[A-Z]\\b"
  while IFS= read -r path; do
    [[ -f "$path" ]] || continue
    case "$path" in
      scripts/pre_push_research_safety_check.sh|twin_engine/privacy.py)
        continue
        ;;
    esac
    if file "$path" | grep -qiE 'text|json|python|html|xml|yaml|script|markdown'; then
      if grep -InE "$phi_regex" "$path" >/tmp/research_safety_grep.$$ 2>/dev/null; then
        staged_phi+="$(cat /tmp/research_safety_grep.$$)"$'\n'
      fi
      rm -f /tmp/research_safety_grep.$$
    fi
  done <<< "$staged_files"
  if [[ -n "$staged_phi" ]]; then
    echo "ERROR: staged text contains possible direct identifiers or PHI markers:" >&2
    echo "$staged_phi" >&2
    exit 1
  fi
fi

echo "Running Django checks and tests..."
"$python_bin" manage.py check
"$python_bin" manage.py test mmportal.tests.test_settings_hosts
"$python_bin" manage.py test twin_engine.tests
"$python_bin" manage.py test clinic.tests.test_patient_crud clinic.tests.test_security_and_ux simulator.tests.test_twin_pr1_integration simulator.tests.test_simulation

echo "Research safety check passed."