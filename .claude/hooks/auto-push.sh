#!/bin/bash
# Auto-commit + push after Edit/Write/NotebookEdit inside this repo.
# Never blocks the agent's flow: all failure paths print a warning and exit 0.
#
# Fixed 2026-10-11:
#   * Paths are derived from this script's own location instead of being hard-coded.
#     The repo moved to ~/Desktop/校外/ISEF, so the old REPO="/Users/ryanwang/Desktop/ISEF"
#     made `cd` fail and the hook exit 0 silently - both agents stopped pushing.
#   * Pushes the *current branch's* upstream instead of a hard-coded `main`.
#     On crosscheck-and-buy-list-v1 the old script committed to the feature branch
#     but pushed `main`, so nothing ever reached GitHub.
#
# Set AUTO_PUSH_DRY_RUN=1 to print what would happen without committing or pushing:
#   echo '{"tool_input":{"file_path":"<some file in this repo>"}}' | AUTO_PUSH_DRY_RUN=1 bash .claude/hooks/auto-push.sh

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOG="$(dirname "${BASH_SOURCE[0]}")/auto-push.log"

input="$(cat)"
file_path="$(printf '%s' "$input" | jq -r '.tool_input.file_path // .tool_input.notebook_path // empty' 2>/dev/null)"

# Only act on edits inside this repo.
case "$file_path" in
  "$REPO"/*) ;;
  *) exit 0 ;;
esac

cd "$REPO" || exit 0

warn() {
  echo "{\"systemMessage\": \"auto-push: $1\"}"
  echo "$(date '+%Y-%m-%d %H:%M:%S') $1" >> "$LOG"
  exit 0
}

branch="$(git rev-parse --abbrev-ref HEAD 2>/dev/null)"
[ -n "$branch" ] || warn "cannot determine current branch, skipped auto-push"
[ "$branch" = "HEAD" ] && warn "detached HEAD, skipped auto-push"

upstream="$(git rev-parse --abbrev-ref --symbolic-full-name '@{u}' 2>/dev/null)"
if [ -n "$upstream" ]; then
  remote="${upstream%%/*}"
else
  remote="origin"
fi

if [ "${AUTO_PUSH_DRY_RUN:-0}" = "1" ]; then
  echo "auto-push dry run"
  echo "  repo     : $REPO"
  echo "  branch   : $branch"
  echo "  upstream : ${upstream:-<none>}"
  echo "  remote   : $remote"
  echo "  changes  :"
  git status --porcelain | sed 's/^/    /'
  exit 0
fi

if [ -z "$(git status --porcelain)" ]; then
  exit 0
fi

if [ -n "$upstream" ]; then
  if ! git fetch "$remote" >>"$LOG" 2>&1; then
    warn "git fetch failed, skipped auto-push (changes stay local/uncommitted)"
  fi

  if ! git pull --rebase >>"$LOG" 2>&1; then
    git rebase --abort >>"$LOG" 2>&1
    warn "rebase against $upstream failed (conflict), skipped auto-push - resolve manually"
  fi
fi

changed="$(git status --porcelain | awk '{print $2}' | xargs -n1 basename 2>/dev/null | paste -sd, -)"
timestamp="$(date '+%Y-%m-%d %H:%M:%S')"

git add -A >>"$LOG" 2>&1

if ! git commit -m "Auto-commit: ${changed:-changes} (${timestamp})" >>"$LOG" 2>&1; then
  warn "git commit failed, changes left staged - check manually"
fi

if ! git push ${upstream:+-u} "$remote" "$branch" >>"$LOG" 2>&1; then
  # Remote moved again between rebase and push - retry once.
  if git pull --rebase >>"$LOG" 2>&1 && git push "$remote" "$branch" >>"$LOG" 2>&1; then
    exit 0
  fi
  warn "git push failed after retry - commit is local only, push manually"
fi

exit 0
