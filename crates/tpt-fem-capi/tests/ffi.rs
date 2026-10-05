//! Exercises the `extern "C"` surface the way a C caller would.

use std::ffi::{c_char, c_void, CStr, CString};
use std::ptr;

use tpt_fem_capi::*;

fn last_error() -> String {
    unsafe { CStr::from_ptr(tpt_last_error_message()) }
        .to_string_lossy()
        .into_owned()
}

fn unit_box(n: usize) -> *mut TptMesh {
    let (min, max, n) = ([0.0; 3], [1.0; 3], [n; 3]);
    let mut mesh = ptr::null_mut();
    assert_eq!(
        unsafe { tpt_mesh_box(min.as_ptr(), max.as_ptr(), n.as_ptr(), &mut mesh) },
        TPT_OK
    );
    assert!(!mesh.is_null());
    mesh
}

/// Node ids on the plane `axis == coord`.
fn plane(mesh: *const TptMesh, axis: usize, coord: f64) -> Vec<usize> {
    let mut count = 0usize;
    unsafe {
        assert_eq!(
            tpt_mesh_nodes_on_plane(mesh, axis, coord, 1e-9, ptr::null_mut(), 0, &mut count),
            TPT_OK
        );
        let mut buf = vec![0usize; count];
        assert_eq!(
            tpt_mesh_nodes_on_plane(mesh, axis, coord, 1e-9, buf.as_mut_ptr(), count, &mut count),
            TPT_OK
        );
        buf
    }
}

fn node_x(mesh: *const TptMesh, i: usize) -> f64 {
    let mut c = [0.0; 3];
    assert_eq!(unsafe { tpt_mesh_coords(mesh, i, c.as_mut_ptr()) }, TPT_OK);
    c[0]
}

#[test]
fn poisson_linear_solution() {
    let mesh = unit_box(3);
    let mut bcs: Vec<TptNodeBc> = plane(mesh, 0, 0.0)
        .into_iter()
        .map(|node| TptNodeBc { node, value: 0.0 })
        .collect();
    bcs.extend(plane(mesh, 0, 1.0).into_iter().map(|node| TptNodeBc { node, value: 1.0 }));

    let mut field = ptr::null_mut();
    let rc = unsafe {
        tpt_solve_poisson(mesh, 1.0, 2, 0.0, None, ptr::null_mut(), bcs.as_ptr(), bcs.len(), &mut field)
    };
    assert_eq!(rc, TPT_OK, "{}", last_error());
    unsafe {
        let n = tpt_mesh_node_count(mesh);
        assert_eq!(tpt_field_len(field), n);
        assert_eq!(tpt_field_dim(field), 1);
        let u = std::slice::from_raw_parts(tpt_field_data(field), n);
        for (i, ui) in u.iter().enumerate() {
            assert!((ui - node_x(mesh, i)).abs() < 1e-8, "node {i}: {ui}");
        }
        tpt_field_free(field);
        tpt_mesh_free(mesh);
    }
}

extern "C" fn counting_source(_: f64, _: f64, _: f64, user: *mut c_void) -> f64 {
    unsafe { *(user as *mut usize) += 1 };
    1.0
}

#[test]
fn poisson_callback_source_is_called() {
    let mesh = unit_box(2);
    let bcs: Vec<TptNodeBc> = plane(mesh, 0, 0.0)
        .into_iter()
        .map(|node| TptNodeBc { node, value: 0.0 })
        .collect();
    let mut calls = 0usize;
    let mut field = ptr::null_mut();
    let rc = unsafe {
        tpt_solve_poisson(
            mesh, 1.0, 2, 0.0, Some(counting_source),
            &mut calls as *mut usize as *mut c_void,
            bcs.as_ptr(), bcs.len(), &mut field,
        )
    };
    assert_eq!(rc, TPT_OK, "{}", last_error());
    assert!(calls > 0);
    unsafe {
        tpt_field_free(field);
        tpt_mesh_free(mesh);
    }
}

#[test]
fn elasticity_and_modal_on_box() {
    let mesh = unit_box(2);
    let model = CString::new("3d").unwrap();
    let bcs: Vec<TptDofBc> = plane(mesh, 0, 0.0)
        .into_iter()
        .flat_map(|node| (0..3).map(move |component| TptDofBc { node, component, value: 0.0 }))
        .collect();
    unsafe {
        let mut field = ptr::null_mut();
        let rc = tpt_solve_elasticity(mesh, model.as_ptr(), 200e9, 0.3, 2, bcs.as_ptr(), bcs.len(), &mut field);
        assert_eq!(rc, TPT_OK, "{}", last_error());
        assert_eq!(tpt_field_dim(field), 3);
        assert_eq!(tpt_field_len(field), 3 * tpt_mesh_node_count(mesh));
        tpt_field_free(field);

        let mut modal = ptr::null_mut();
        let rc = tpt_solve_modal(mesh, model.as_ptr(), 200e9, 0.3, 7800.0, 2, 2, bcs.as_ptr(), bcs.len(), &mut modal);
        assert_eq!(rc, TPT_OK, "{}", last_error());
        assert_eq!(tpt_modal_count(modal), 2);
        assert!(tpt_modal_omega2(modal, 0) > 0.0);
        assert!(tpt_modal_omega2(modal, 99).is_nan());
        assert_eq!(tpt_modal_shape_len(modal, 0), 3 * tpt_mesh_node_count(mesh));
        assert!(tpt_modal_shape(modal, 99).is_null());
        tpt_modal_free(modal);
        tpt_mesh_free(mesh);
    }
}

#[test]
fn topopt_small_cantilever() {
    let mut t = ptr::null_mut();
    let rc = unsafe { tpt_topopt_cantilever(6, 3, 0.5, 3.0, 1.5, 5, &mut t) };
    assert_eq!(rc, TPT_OK, "{}", last_error());
    unsafe {
        assert_eq!((tpt_topopt_nx(t), tpt_topopt_ny(t)), (6, 3));
        let d = std::slice::from_raw_parts(tpt_topopt_densities(t), 18);
        assert!(d.iter().all(|v| (0.0..=1.0 + 1e-9).contains(v)));
        assert!(tpt_topopt_compliance_len(t) > 0);
        tpt_topopt_free(t);
    }
}

#[test]
fn invalid_arguments_report_errors() {
    unsafe {
        let mut mesh = ptr::null_mut();
        assert_eq!(tpt_mesh_box(ptr::null(), ptr::null(), ptr::null(), &mut mesh), TPT_ERR_INVALID_ARG);
        assert!(last_error().contains("NULL"));

        let path = CString::new("/definitely/not/here.msh").unwrap();
        assert_eq!(tpt_mesh_load(path.as_ptr(), &mut mesh), TPT_ERR_IO);
        assert!(mesh.is_null());

        let mesh = unit_box(1);
        let bad = CString::new("sphere").unwrap();
        let mut f = ptr::null_mut();
        assert_eq!(
            tpt_solve_elasticity(mesh, bad.as_ptr(), 1.0, 0.3, 2, ptr::null(), 0, &mut f),
            TPT_ERR_INVALID_ARG
        );
        assert!(last_error().contains("unknown elasticity model"));
        assert_eq!(tpt_mesh_coords(mesh, 10_000, [0.0; 3].as_mut_ptr()), TPT_ERR_INVALID_ARG);

        let mut t = ptr::null_mut();
        assert_eq!(tpt_topopt_cantilever(0, 3, 0.5, 3.0, 1.5, 5, &mut t), TPT_ERR_INVALID_ARG);

        // NULL frees and accessors are no-ops.
        tpt_mesh_free(ptr::null_mut());
        tpt_field_free(ptr::null_mut());
        assert_eq!(tpt_field_len(ptr::null()), 0);
        tpt_mesh_free(mesh);
    }
}

#[test]
fn version_is_exposed() {
    let v = unsafe { CStr::from_ptr(tpt_version() as *const c_char) };
    assert_eq!(v.to_str().unwrap(), env!("CARGO_PKG_VERSION"));
}
