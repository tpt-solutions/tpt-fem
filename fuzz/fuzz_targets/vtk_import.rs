#![no_main]
//! Fuzz the promoted `tpt-fem-io-vtk` VTK reader (`read_vtk`) against arbitrary
//! bytes. The reader is the untrusted-input surface for `.vtk`/`.vtu` import, so
//! it must not panic on malformed input — it should return `Err` instead.
use libfuzzer_sys::fuzz_target;
use std::io::Write;

fuzz_target!(|data: &[u8]| {
    let dir = std::env::temp_dir();
    let path = dir.join(format!("tpt_fem_fuzz_vtk_{}.vtk", std::process::id()));
    if let Ok(mut f) = std::fs::File::create(&path) {
        let _ = f.write_all(data);
        // The reader must either succeed or return an error — never panic.
        let _ = tpt_fem_io_vtk::read_vtk(&path);
    }
    let _ = std::fs::remove_file(&path);
});
