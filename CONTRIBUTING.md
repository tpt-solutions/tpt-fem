# Contributing to tpt-fem

Thanks for your interest in improving the tpt-fem workspace. This document
covers how to report issues and what we expect from a good bug report or
feature request.

## Reporting issues

Before opening a new issue, please search the existing issues to avoid
duplicates.

### Bug reports

A good bug report includes:

- A clear, descriptive title.
- The exact steps to reproduce the problem.
- The expected behavior and what happened instead.
- The environment: OS, Rust toolchain version (`rustc --version`,
  `cargo --version`), and which crate(s) are affected.
- A minimal reproducible example where possible (a small test or snippet).
- Any relevant compiler output, logs, or panics, ideally trimmed to the
  relevant lines.

### Feature requests

A good feature request includes:

- The problem you are trying to solve and your use case.
- A description of the proposed change or new functionality.
- Whether it fits within an existing crate's responsibilities, or whether a
  new crate would be needed (see the crate dependency DAG in the repository
  for layering constraints).
- Any alternatives you considered.

## Issue labels and triage

Issues are triaged using labels for area (e.g. `crate: thermal`, `crate:
sparse`), kind (`bug`, `enhancement`, `docs`, `question`), and priority. Please
do not apply labels yourself unless you are a maintainer; set the issue body
and let maintainers classify it.

## Scope and layering

The workspace is a layered DAG. When proposing a change, keep in mind that lower
layers never depend on higher layers and cross-layer cycles are not allowed. If
your issue touches the crate structure, mention which crates are involved.

## License

Contributions are accepted under the same dual `MIT OR Apache-2.0` license as
the project.
