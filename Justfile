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
