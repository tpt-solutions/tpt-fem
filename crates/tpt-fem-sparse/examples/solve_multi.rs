use tpt_fem_sparse::Coo;

fn main() {
    let mut c = Coo::new();
    c.push(0, 0, 2.0);
    c.push(0, 1, 1.0);
    c.push(1, 0, 1.0);
    c.push(1, 1, 3.0);

    // Solve against two right-hand sides with a single (implicit) factorization.
    let sols = tpt_fem_sparse::solve_multi(&c, &[vec![3.0, 5.0], vec![1.0, 1.0]]).expect("solve");

    // First RHS [3, 5] -> [0.8, 1.4]; second RHS [1, 1] -> [0.4, 0.2].
    assert!((sols[0][0] - 0.8).abs() < 1e-10);
    assert!((sols[0][1] - 1.4).abs() < 1e-10);
    assert!((sols[1][0] - 0.4).abs() < 1e-10);
    assert!((sols[1][1] - 0.2).abs() < 1e-10);
    println!("Solutions = {:?}", sols);
}
