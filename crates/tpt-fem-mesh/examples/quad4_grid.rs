use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let nx = 2usize;
    let ny = 2usize;
    let mut b = MeshBuilder::new();
    let mut ids = vec![vec![0usize; ny + 1]; nx + 1];
    for i in 0..=nx {
        for j in 0..=ny {
            let x = i as f64 / nx as f64;
            let y = j as f64 / ny as f64;
            ids[i][j] = b.add_node(vec![x, y]);
        }
    }
    for i in 0..nx {
        for j in 0..ny {
            b.add_element(
                CellType::Quad,
                vec![ids[i][j], ids[i + 1][j], ids[i + 1][j + 1], ids[i][j + 1]],
            );
        }
    }
    let mesh = b.build();

    println!(
        "node_count   = {} (expected {})",
        mesh.node_count(),
        (nx + 1) * (ny + 1)
    );
    println!(
        "element_count = {} (expected {})",
        mesh.element_count(),
        nx * ny
    );

    assert_eq!(mesh.node_count(), (nx + 1) * (ny + 1));
    assert_eq!(mesh.element_count(), nx * ny);
    assert_eq!(mesh.node_coords(ids[0][0]), &[0.0, 0.0]);
    assert_eq!(mesh.node_coords(ids[nx][ny]), &[1.0, 1.0]);
}
