//! Regenerates `include/tpt_fem.h` with cbindgen. A failure only warns, so a
//! cbindgen hiccup never breaks the library build (CI diffs the header).

fn main() {
    let dir = std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR is set");
    println!("cargo:rerun-if-changed=src");
    println!("cargo:rerun-if-changed=cbindgen.toml");
    let config = cbindgen::Config::from_file(format!("{dir}/cbindgen.toml"))
        .expect("cbindgen.toml is valid");
    match cbindgen::Builder::new()
        .with_crate(&dir)
        .with_config(config)
        .generate()
    {
        Ok(bindings) => {
            bindings.write_to_file(format!("{dir}/include/tpt_fem.h"));
        }
        Err(e) => println!("cargo:warning=cbindgen failed: {e}"),
    }
}
