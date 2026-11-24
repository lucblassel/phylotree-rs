#!/usr/bin/env bash
# Run the same checks as .github/workflows/rustcheck.yaml locally.
#
# Usage:
#   scripts/ci.sh              run all steps
#   scripts/ci.sh fmt clippy   run only the named steps (fmt check test clippy doc)
set -euo pipefail

cd "$(dirname "$0")/.."

steps=("$@")
if [ ${#steps[@]} -eq 0 ]; then
  steps=(fmt check test clippy doc)
fi

for step in "${steps[@]}"; do
  case "$step" in
    fmt)
      echo "== cargo fmt --all -- --check"
      cargo fmt --all -- --check
      ;;
    check)
      # Benches are excluded: shared runners give non-reproducible timings.
      # Run them locally with `cargo bench` when needed.
      echo "== cargo check --lib --bins --tests"
      cargo check --lib --bins --tests
      echo "== cargo check --no-default-features"
      cargo check --no-default-features
      # All features together prove the feature set is mutually compatible.
      echo "== cargo check --all-features"
      cargo check --all-features
      ;;
    test)
      echo "== cargo test --lib --bins --tests"
      cargo test --lib --bins --tests
      echo "== cargo test --doc"
      cargo test --doc
      # Long-running deep-tree tests are #[ignore]d by default; run them here.
      echo "== cargo test --lib --tests --bins -- --ignored deep_"
      cargo test --lib --tests --bins -- --ignored deep_
      ;;
    clippy)
      echo "== cargo clippy --lib --bins --tests -- -D warnings"
      cargo clippy --lib --bins --tests -- -D warnings
      ;;
    doc)
      echo "== cargo doc --no-deps (RUSTDOCFLAGS=-D warnings)"
      RUSTDOCFLAGS="-D warnings" cargo doc --no-deps
      ;;
    *)
      echo "unknown step: $step (expected one of: fmt check test clippy doc)" >&2
      exit 1
      ;;
  esac
done
