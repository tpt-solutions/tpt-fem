#![no_main]
libfuzzer_sys::fuzz_target!(|data: &[u8]| {
    // Abaqus `.inp` import must reject malformed decks via `Result`, not panic.
    if let Ok(s) = std::str::from_utf8(data) {
        let _ = tpt_fem_io_abaqus::read_inp(s);
    }
});
