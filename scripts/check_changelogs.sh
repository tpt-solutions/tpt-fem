#!/usr/bin/env bash
# Fail if a crate's `src/` changed relative to a base ref without a matching
# CHANGELOG.md update. Usage: scripts/check_changelogs.sh [base-ref]
# (default: origin/master). Used by CI on pull requests.
set -euo pipefail
base="${1:-origin/master}"
changed="$(git diff --name-only "$base"...HEAD)"
status=0
for dir in crates/*/; do
  crate="$(basename "$dir")"
  [ -f "$dir/CHANGELOG.md" ] || continue
  if grep -q "^crates/$crate/src/" <<<"$changed" \
     && ! grep -q "^crates/$crate/CHANGELOG.md$" <<<"$changed"; then
    echo "MISSING CHANGELOG entry: $crate (src/ changed, CHANGELOG.md did not)"
    status=1
  fi
done
[ "$status" -eq 0 ] && echo "changelogs are up to date for all changed crates"
exit "$status"
