# tpt-fem-wasm

In-browser SIMP topology-optimization demo: the `tpt-fem-topopt` solver compiled
to WebAssembly with a small static page (`www/`). Dev-only (`publish = false`,
excluded from the Cargo workspace like `tpt-fem-py`).

```sh
wasm-pack build --target web --out-dir www/pkg
python -m http.server -d www 8080   # then open http://localhost:8080
```

The solve runs on the main thread using the in-house dense LU, so the grid is
capped at 40 × 16 elements (the default 30 × 10 takes a couple of seconds; the cap takes tens of seconds). `www/` is static and can be hosted as-is (e.g.
GitHub Pages) after building `www/pkg`.
