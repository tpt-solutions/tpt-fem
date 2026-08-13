#![no_main]
libfuzzer_sys::fuzz_target!(|data: &[u8]| {
    // Gmsh `.msh` import must reject malformed input via `Result`, not panic.
    let _ = tpt_fem_mesh::Mesh::from_msh_bytes(data);
});
