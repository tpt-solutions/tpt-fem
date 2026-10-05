# Justfile for the tpt-fem workspace.
# Run `just` to see the available recipes, or `just <recipe>`.

# List available recipes.
default:
    @just --list

# Build the whole workspace with every feature enabled.
build:
    cargo build --workspace --all-features

# Run the full test suite with every feature enabled.
test:
    cargo test --workspace --all-features

# Check formatting without modifying files.
fmt:
    cargo fmt --all --check

# Format all sources in place.
fmt-fix:
    cargo fmt --all

# Lint with clippy, treating warnings as errors.
clippy:
    cargo clippy --workspace --all-targets --all-features -- -D warnings

# Run the license/dependency audit.
deny:
    cargo deny check

# Comprehensive verification gate used in CI.
verify: fmt clippy deny test

# Run the prelude-only thermal solve example.
example-thermal:
    cargo run -p tpt-fem --example thermal_solve

# Run the ignored (heavy) 3-D MMS convergence test.
mms-3d:
    cargo test -p tpt-fem-thermal --test mms_convergence -- --ignored

# Build the command-line driver.
cli:
    cargo build -p tpt-fem-cli

# Verify the excluded crates (py/capi/wasm) pin the workspace's tpt-fem* versions.
pins:
    python3 scripts/check_pins.py

# Verify every crate whose src/ changed vs a base ref also updated its CHANGELOG.
changelogs base="origin/master":
    scripts/check_changelogs.sh {{base}}

# Check the declared MSRV (keep in sync with rust-version in Cargo.toml).
msrv:
    cargo +1.85.0 check --workspace --all-targets

# Validate every publishable crate's package contents.
package:
    for d in crates/tpt-fem*/; do n=$(basename $d); case $n in tpt-fem-py|tpt-fem-capi|tpt-fem-wasm) continue;; esac; cargo package -p $n --list --allow-dirty > /dev/null || exit 1; done

# Every gate in PUBLISHING.md that can run locally, before the first `cargo publish`.
release-check: verify pins msrv package
    RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
    cargo test -p tpt-fem --test manifest_drift
    cargo build --manifest-path crates/tpt-fem-capi/Cargo.toml
    cargo check --manifest-path crates/tpt-fem-py/Cargo.toml
    cargo check --manifest-path crates/tpt-fem-wasm/Cargo.toml --target wasm32-unknown-unknown
