#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

python_bin="${PYTHON_BIN:-$repo_root/.venv/bin/python}"
if [[ ! -x "$python_bin" && -x "$repo_root/venv/bin/python" ]]; then
  python_bin="$repo_root/venv/bin/python"
fi
[[ -x "$python_bin" ]] || fail "Python interpreter not found or not executable: $python_bin"

base_ref="${BMYCURE4MM_SAFETY_BASE_REF:-}"
head_ref="${BMYCURE4MM_SAFETY_HEAD_REF:-}"
safety_mode="local"
resolved_base=""
resolved_head=""

if [[ -n "$base_ref" || -n "$head_ref" ]]; then
  [[ -n "$base_ref" && -n "$head_ref" ]] || \
    fail "CI safety mode requires both BMYCURE4MM_SAFETY_BASE_REF and BMYCURE4MM_SAFETY_HEAD_REF."
  [[ "$base_ref" != -* && "$head_ref" != -* ]] || fail "Safety refs must not begin with '-'."
  resolved_base="$(git rev-parse --verify "${base_ref}^{commit}" 2>/dev/null)" || \
    fail "Invalid safety base ref."
  resolved_head="$(git rev-parse --verify "${head_ref}^{commit}" 2>/dev/null)" || \
    fail "Invalid safety head ref."
  safety_mode="ci"
fi

echo "== Research safety check ($safety_mode mode) =="

required_ignored=("local_private/.safety-probe" "media/.safety-probe" "db.sqlite3" ".env")
for path in "${required_ignored[@]}"; do
  git check-ignore -q "$path" || fail "Expected ignored path is not ignored: $path"
done

tracked_sensitive="$({
  git ls-files | grep -E \
    '(^|/)(local_private|media|logs)(/|$)|(^|/)\.env($|\.)|(^|/)db\.sqlite3$|\.(sqlite3|db|log|pem|key|p12|pdf)$' || true
} | grep -v -E '^\.env\.example$|^modules/binding_visualizer/sources/5LF3_structure_report\.pdf$' || true)"
if [[ -n "$tracked_sensitive" ]]; then
  echo "ERROR: tracked sensitive or prohibited artifact paths detected:" >&2
  printf '%s\n' "$tracked_sensitive" >&2
  exit 1
fi

changed_paths=()
if [[ "$safety_mode" == "ci" ]]; then
  while IFS= read -r -d '' path; do
    changed_paths+=("$path")
  done < <(git diff --no-ext-diff --name-only -z --diff-filter=ACMRT "$resolved_base" "$resolved_head" --)
else
  while IFS= read -r -d '' path; do
    changed_paths+=("$path")
  done < <(git diff --cached --no-ext-diff --name-only -z --diff-filter=ACMRT --)
fi

direct_name_one="Ros""sana"
direct_name_two="Agu""eci"
tax_marker="Codice"" Fiscale"
birth_marker="Data"" di nascita"
phi_regex="\\b(${direct_name_one}|${direct_name_two})\\b|${tax_marker}|${birth_marker}|\\b[A-Z]{6}[0-9]{2}[A-Z][0-9]{2}[A-Z][0-9]{3}[A-Z]\\b"
private_key_marker="-----BEGIN [A-Z ]*PRIVATE"" KEY-----"
github_token_marker="gh""[pousr]_[A-Za-z0-9]{20,}"
aws_access_marker="AK""IA[0-9A-Z]{16}"
credential_regex="${private_key_marker}|${github_token_marker}|${aws_access_marker}"
absolute_path_regex='/Users/[A-Za-z0-9._-]+/|/Volumes/[A-Za-z0-9._-]+/|/home/[A-Za-z0-9._-]+/|[A-Za-z]:\\Users\\'

path_matches() {
  local path="$1"
  local regex="$2"
  if [[ "$safety_mode" == "ci" ]]; then
    git grep -Eq -e "$regex" "$resolved_head" -- "$path" >/dev/null 2>&1
  else
    git grep --cached -Eq -e "$regex" -- "$path" >/dev/null 2>&1
  fi
}

scan_blob_for_findings() {
  local path="$1"
  if path_matches "$path" "$credential_regex"; then
    echo "ERROR: possible credential or private key detected in tracked content: $path" >&2
    return 1
  fi
  case "$path" in
    scripts/pre_push_research_safety_check.sh|twin_engine/privacy.py)
      ;;
    *)
      if path_matches "$path" "$phi_regex"; then
        echo "ERROR: possible direct identifier or PHI marker detected in tracked content: $path" >&2
        return 1
      fi
      ;;
  esac
  case "$path" in
    docs/*|*/tests/*|tests/*|scripts/pre_push_research_safety_check.sh)
      ;;
    *.py|*.sh|*.yml|*.yaml|Dockerfile|docker-compose.yml)
      if path_matches "$path" "$absolute_path_regex"; then
        echo "ERROR: absolute local filesystem path detected in runtime code: $path" >&2
        return 1
      fi
      ;;
  esac
}

if [[ ${#changed_paths[@]} -gt 0 ]]; then
  for path in "${changed_paths[@]}"; do
    if printf '%s\n' "$path" | grep -Eq \
      '(^|/)(local_private|media|logs)(/|$)|(^|/)\.env($|\.)|(^|/)db\.sqlite3$|\.(sqlite3|db|log|pem|key|p12|pdf|png|jpg|jpeg|tif|tiff)$'; then
      fail "Changed sensitive or binary artifact path detected: $path"
    fi
    scan_blob_for_findings "$path"
  done
fi

# Repository-wide high-confidence content invariants query Git's tracked index;
# ignored private directories are never traversed and matched content is hidden.
credential_files="$(git grep --cached -IlE -e "$credential_regex" -- . 2>/dev/null || true)"
if [[ -n "$credential_files" ]]; then
  fail "Possible credential or private key detected in tracked content: $credential_files"
fi
phi_files="$(git grep --cached -IlE -e "$phi_regex" -- . 2>/dev/null | \
  grep -v -E '^scripts/pre_push_research_safety_check\.sh$|^twin_engine/privacy\.py$' || true)"
if [[ -n "$phi_files" ]]; then
  fail "Possible direct identifier or PHI marker detected in tracked content: $phi_files"
fi

compile_root="$(mktemp -d "${TMPDIR:-/tmp}/bmycure4mm-compile.XXXXXX")"
cleanup_compile_root() {
  case "$compile_root" in
    "${TMPDIR:-/tmp}"/bmycure4mm-compile.*) rm -rf -- "$compile_root" ;;
    *) fail "Refusing to remove unexpected compile directory." ;;
  esac
}
trap cleanup_compile_root EXIT
export PYTHONPYCACHEPREFIX="$compile_root"
while IFS= read -r -d '' path; do
  "$python_bin" -m py_compile "$path"
done < <(git ls-files -z '*.py')

if [[ "${BMYCURE4MM_SAFETY_SCAN_ONLY:-0}" == "1" ]]; then
  echo "Research safety scan passed (test scan-only mode)."
  exit 0
fi

echo "Running Django checks and focused safety tests..."
"$python_bin" manage.py check
"$python_bin" manage.py test mmportal.tests.test_settings_hosts
"$python_bin" manage.py test twin_engine.tests
"$python_bin" manage.py test \
  clinic.tests.test_patient_crud \
  clinic.tests.test_security_and_ux \
  simulator.tests.test_twin_pr1_integration \
  simulator.tests.test_simulation

echo "Research safety check passed."
