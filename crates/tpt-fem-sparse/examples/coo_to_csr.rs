use tpt_fem_sparse::Coo;

fn main() {
    // Assemble [[2, 1], [1, 3]] by writing each entry twice, then summing.
    let mut c = Coo::new();
    c.push(0, 0, 1.0);
    c.push(0, 0, 1.0);
    c.push(0, 1, 1.0);
    c.push(1, 0, 1.0);
    c.push(1, 1, 1.5);
    c.push(1, 1, 1.5);
    let csr = c.to_csr();

    assert_eq!(csr.nnz(), 4);
    assert_eq!(csr.row_ptrs, vec![0, 2, 4]);
    assert_eq!(csr.col_ind, vec![0, 1, 0, 1]);
    assert_eq!(csr.values, vec![2.0, 1.0, 1.0, 3.0]);
    println!("CSR nnz = {}, values = {:?}", csr.nnz(), csr.values);
}
