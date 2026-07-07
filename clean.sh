#!/usr/bin/env bash
set -euo pipefail

usage() {
    printf 'Usage: %s [--dry-run] [--tracked-generated]\n' "$0"
    printf '\n'
    printf 'Remove ignored generated files from the working tree.\n'
    printf '\n'
    printf 'Options:\n'
    printf '  --dry-run             Show what would be removed or restored.\n'
    printf '  --tracked-generated   Restore tracked generated cache/build files to HEAD.\n'
    printf '  -h, --help            Show this help text.\n'
}

dry_run=0
tracked_generated=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            dry_run=1
            ;;
        --tracked-generated)
            tracked_generated=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'Error: unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
    shift
done

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    printf 'Error: not inside a git repository.\n' >&2
    exit 1
fi

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"

pathspecs=(
    ':(glob)**/__pycache__/**'
    ':(glob)**/*.py[cod]'
    ':(glob)**/*$py.class'
    'build/**'
    'dist/**'
    '*.egg-info/**'
    '.pytest_cache/**'
    '.coverage'
    'htmlcov/**'
)

if [[ "$dry_run" -eq 1 ]]; then
    printf 'Ignored/generated files that would be removed:\n'
    git clean -ndX -- "${pathspecs[@]}"
else
    git clean -fdX -- "${pathspecs[@]}"
fi

if [[ "$tracked_generated" -eq 1 ]]; then
    if [[ "$dry_run" -eq 1 ]]; then
        printf '\nTracked generated files that would be restored:\n'
        git diff --name-only -- "${pathspecs[@]}"
    else
        git restore -- "${pathspecs[@]}"
    fi
fi
