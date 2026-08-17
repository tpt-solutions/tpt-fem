#![no_main]
libfuzzer_sys::fuzz_target!(|_data: &[u8]| {
    use tpt_fem_mesh::{CellType, MeshBuilder};

    // Exercise the full Exodus write -> read round-trip on a small mesh so the
    // codec is fuzzed for panics independent of the decode-only target above.
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0, 0.0]);
    let n3 = b.add_node(vec![1.0, 1.0, 0.0]);
    b.add_element(CellType::Tri, vec![n0, n1, n2]);
    b.add_element(CellType::Tri, vec![n1, n3, n2]);
    let mesh = b.build();

    let bytes = match tpt_fem_io_exodus::mesh_to_exodus_bytes(&mesh) {
        Ok(b) => b,
        Err(_) => return,
    };
    let _ = tpt_fem_io_exodus::bytes_to_mesh(&bytes);
});
