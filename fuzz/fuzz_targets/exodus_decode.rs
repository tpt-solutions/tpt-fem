#![no_main]
libfuzzer_sys::fuzz_target!(|data: &[u8]| {
    // The hand-rolled NetCDF-3 codec must never panic on untrusted bytes.
    let _ = tpt_fem_io_exodus::bytes_to_mesh(data);
});
