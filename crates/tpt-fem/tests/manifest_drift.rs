// d:\Programming\1PRODUCTION\Open Source\tpt-fem\crates\tpt-fem\tests\manifest_drift.rs
//! Manifest drift guards (todo.md Phase 14c).
//!
//! These tests mechanically catch the class of manifest/merge bug that hit
//! `tpt-fem-modal` in Phase 14a: a crate directory that is a workspace member
//! but missing from `[workspace.dependencies]` (or vice versa). They parse the
//! root `Cargo.toml` directly, so they fail in seconds — long before a full
//! workspace build would notice.

use std::collections::BTreeSet;
use std::path::PathBuf;

fn root_manifest() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("Cargo.toml")
}

/// The `[workspace] members` list from the root manifest.
fn workspace_members(root: &str) -> Vec<String> {
    let value: toml::Value = toml::from_str(root).expect("root Cargo.toml parses");
    value["workspace"]["members"]
        .as_array()
        .expect("workspace.members is an array")
        .iter()
        .map(|v| v.as_str().expect("member paths are strings").to_string())
        .collect()
}

/// The keys of the `[workspace.dependencies]` table.
fn workspace_dependencies(root: &str) -> BTreeSet<String> {
    let value: toml::Value = toml::from_str(root).expect("root Cargo.toml parses");
    value["workspace"]["dependencies"]
        .as_table()
        .expect("workspace.dependencies is a table")
        .keys()
        .cloned()
        .collect()
}

#[test]
fn every_member_appears_in_workspace_dependencies() {
    let root = std::fs::read_to_string(root_manifest()).expect("root Cargo.toml is readable");
    let members = workspace_members(&root);
    assert!(!members.is_empty(), "workspace has members");
    let deps = workspace_dependencies(&root);
    for m in &members {
        let name = m.rsplit('/').next().unwrap();
        // The umbrella and CLI are members but are not consumed as workspace
        // dependencies by other crates via path lookups... except they are
        // declared there too. Every member must have an entry so internal
        // `{ workspace = true }` references can never dangle.
        assert!(
            deps.contains(name),
            "workspace member '{m}' ({name}) is missing from [workspace.dependencies]"
        );
    }
}

#[test]
fn every_workspace_dependency_is_a_member() {
    let root = std::fs::read_to_string(root_manifest()).expect("root Cargo.toml is readable");
    let members: BTreeSet<String> = workspace_members(&root).into_iter().collect();
    for dep in workspace_dependencies(&root) {
        let dir = format!("crates/{dep}");
        assert!(
            members.contains(&dir),
            "[workspace.dependencies] entry '{dep}' does not correspond to a \
             workspace member ('{dir}' not in workspace.members)"
        );
    }
}

#[test]
fn every_crates_directory_is_a_workspace_member() {
    let root = std::fs::read_to_string(root_manifest()).expect("root Cargo.toml is readable");
    let members: BTreeSet<String> = workspace_members(&root).into_iter().collect();
    let crates_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("crates");
    let mut found = 0;
    for entry in std::fs::read_dir(&crates_dir).expect("crates/ is readable") {
        let entry = entry.expect("crates/ entry readable");
        if !entry.file_type().expect("file type readable").is_dir() {
            continue;
        }
        let name = entry.file_name().to_string_lossy().to_string();
        if !name.starts_with("tpt-fem-") {
            continue;
        }
        // tpt-fem-py is excluded from the workspace (built by maturin).
        if name == "tpt-fem-py" {
            continue;
        }
        found += 1;
        let member = format!("crates/{name}");
        assert!(
            members.contains(&member),
            "crate directory 'crates/{name}' exists on disk but is NOT a \
             workspace member in the root Cargo.toml"
        );
    }
    assert!(found > 0, "found no tpt-fem-* crates under crates/");
}
