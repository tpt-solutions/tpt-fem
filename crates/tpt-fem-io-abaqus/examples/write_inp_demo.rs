use tpt_fem_io_abaqus::write_inp;
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0]);
    b.add_element(CellType::Tri, vec![n0, n1, n2]);
    let mesh = b.build();

    let path = std::env::temp_dir().join("tpt_fem_io_abaqus_demo.inp");
    write_inp(&mesh, &path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    println!("{text}");

    assert!(text.lines().any(|l| l.starts_with("*NODE")));
    assert!(text.contains("*ELEMENT, TYPE=CPS3"));
    // Abaqus uses 1-based ids, so the first node must be written as "1,".
    assert!(text.lines().any(|l| l.starts_with("1,")));
    let _ = std::fs::remove_file(&path);
}
