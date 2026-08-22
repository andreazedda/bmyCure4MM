#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

phase="${1:-all}"
case "$phase" in
  all|repository-hygiene|quality|django-integrity|scientific-regression|dependency-security|container-build)
    ;;
  *)
    echo "ERROR: unknown required-check phase: $phase" >&2
    exit 2
    ;;
esac

if [[ "$phase" != "container-build" ]] && ! command -v uv >/dev/null 2>&1; then
  if [[ -x "$repo_root/venv/bin/uv" ]]; then
    export PATH="$repo_root/venv/bin:$PATH"
  else
    echo "ERROR: uv 0.12.3 is required on PATH." >&2
    exit 1
  fi
fi

initial_status_hash="$(git status --porcelain=v1 -z --untracked-files=all | git hash-object --stdin)"
ci_temp_root="$(mktemp -d "${TMPDIR:-/tmp}/bmycure4mm-required.XXXXXX")"
cleanup_ci_temp_root() {
  case "$ci_temp_root" in
    "${TMPDIR:-/tmp}"/bmycure4mm-required.*) rm -rf -- "$ci_temp_root" ;;
    *)
      echo "ERROR: refusing to remove unexpected CI temporary directory." >&2
      return 1
      ;;
  esac
}
trap cleanup_ci_temp_root EXIT

export BMYCURE4MM_CI_TEMP_ROOT="$ci_temp_root"
export DJANGO_SETTINGS_MODULE="mmportal.settings_ci"
export DJANGO_DEBUG="0"
export DJANGO_SECRET_KEY="ci-only-synthetic-key-not-used-by-settings-module"
export DJANGO_SQLITE_PATH="$ci_temp_root/db/ci.sqlite3"
export DJANGO_MEDIA_ROOT="$ci_temp_root/media"
export DJANGO_LOGS_ROOT="$ci_temp_root/logs"
export CELERY_TASK_ALWAYS_EAGER="1"
export OPENBLAS_NUM_THREADS="1"
export OMP_NUM_THREADS="1"
export MKL_NUM_THREADS="1"
export VECLIB_MAXIMUM_THREADS="1"
export PYTHONHASHSEED="0"
export UV_LINK_MODE="copy"

uv_python_args=()
if [[ -n "${BMYCURE4MM_CI_PYTHON:-}" ]]; then
  uv_python_args=(--python "$BMYCURE4MM_CI_PYTHON")
fi

uv_run() {
  uv run --frozen --extra chemistry "${uv_python_args[@]}" "$@"
}

sync_environment() {
  echo "== locked environment =="
  uv --version
  uv sync --frozen --extra chemistry "${uv_python_args[@]}"
}

run_repository_hygiene() {
  echo "== repository hygiene and research safety =="
  uv_run python -m scripts.check_repository_hygiene
  uv_run python -m unittest tests.test_ci_governance
  PYTHON_BIN="$repo_root/.venv/bin/python" bash scripts/pre_push_research_safety_check.sh
}

run_quality() {
  echo "== formatting, lint, and typing =="
  uv_run ruff format --check .
  uv_run ruff check .
  uv_run mypy
}

run_django_integrity() {
  echo "== Django integrity =="
  uv_run python manage.py check --settings=mmportal.settings_ci
  uv_run python manage.py check --deploy --fail-level WARNING --settings=mmportal.settings_ci_deploy
  uv_run python manage.py makemigrations --check --dry-run --settings=mmportal.settings_ci
  uv_run python manage.py test --settings=mmportal.settings_ci
}

run_scientific_regression() {
  echo "== scientific regression contracts =="
  uv_run python -m scripts.check_numerical_baseline
  uv_run python manage.py test twin_engine simulator clinic --settings=mmportal.settings_ci
}

run_dependency_security() {
  echo "== dependency security =="
  uv_run python -m scripts.audit_dependencies
}

run_container_build() {
  echo "== container build and privacy-safe smoke check =="
  command -v docker >/dev/null 2>&1 || {
    echo "ERROR: docker is required for the container-build phase." >&2
    return 1
  }
  local git_sha image_tag source_label
  git_sha="$(git rev-parse HEAD)"
  image_tag="bmycure4mm:m0-r-ci-${git_sha:0:12}"
  source_label="https://github.com/andreazedda/bmyCure4MM"
  docker build \
    --label "org.opencontainers.image.source=$source_label" \
    --label "org.opencontainers.image.revision=$git_sha" \
    --tag "$image_tag" \
    .
  test "$(docker image inspect --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}' "$image_tag")" = "$git_sha"
  docker run --rm --entrypoint sh "$image_tag" -c '
    test ! -e /app/local_private
    test ! -e /app/db.sqlite3
    test ! -e /app/.env
    test -z "$(find /app/media /app/logs -type f -print -quit)"
  '
  docker run --rm \
    --entrypoint python \
    -e BMYCURE4MM_CI_TEMP_ROOT=/tmp/bmycure4mm-ci \
    "$image_tag" \
    manage.py check --settings=mmportal.settings_ci
}

if [[ "$phase" == "container-build" ]]; then
  run_container_build
else
  sync_environment
  case "$phase" in
    repository-hygiene) run_repository_hygiene ;;
    quality) run_quality ;;
    django-integrity) run_django_integrity ;;
    scientific-regression) run_scientific_regression ;;
    dependency-security) run_dependency_security ;;
    all)
      run_repository_hygiene
      run_quality
      run_django_integrity
      run_scientific_regression
      run_dependency_security
      run_container_build
      ;;
  esac
fi

final_status_hash="$(git status --porcelain=v1 -z --untracked-files=all | git hash-object --stdin)"
if [[ "$initial_status_hash" != "$final_status_hash" ]]; then
  echo "ERROR: required checks modified the Git working tree." >&2
  exit 1
fi

echo "M0_R_LOCAL_REQUIRED_CHECKS = PASS"
