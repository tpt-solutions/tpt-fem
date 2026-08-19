use tpt_fem_sparse::Coo;

fn main() {
    // Solve [[2, 1], [1, 3]] x = [3, 5].
    // Hand solution: x = [0.8, 1.4].
    let mut c = Coo::new();
    c.push(0, 0, 2.0);
    c.push(0, 1, 1.0);
    c.push(1, 0, 1.0);
    c.push(1, 1, 3.0);

    let x = tpt_fem_sparse::solve(&c, &[3.0, 5.0]).expect("solve");
    assert!((x[0] - 0.8).abs() < 1e-10);
    assert!((x[1] - 1.4).abs() < 1e-10);
    println!("Solution x = {:?}", x);
}
