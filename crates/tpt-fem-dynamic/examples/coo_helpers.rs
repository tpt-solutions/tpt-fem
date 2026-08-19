//! Standalone demo of the `coo_scale`, `coo_add`, and `coo_matvec` helpers.
//!
//! Run with: `cargo run -p tpt-fem-dynamic --example coo_helpers`
//!
//! Builds `A = [[2, 1], [1, 3]]` and verifies the sparse helpers against exact
//! dense arithmetic: `2·A == A + A` and `A·[1, 2]ᵀ == [4, 7]ᵀ`.

use tpt_fem_dynamic::{coo_add, coo_matvec, coo_scale};
use tpt_fem_sparse::Coo;

fn main() {
    let mut a = Coo::new();
    a.push(0, 0, 2.0);
    a.push(0, 1, 1.0);
    a.push(1, 0, 1.0);
    a.push(1, 1, 3.0);

    // coo_scale: 2·A = [[4, 2], [2, 6]].
    let two_a = coo_scale(&a, 2.0);
    assert_eq!(two_a.vals, vec![4.0, 2.0, 2.0, 6.0]);
    println!("coo_scale(A, 2.0).vals = {:?}", two_a.vals);

    // coo_add: A + A == 2·A. The raw triplet list is concatenated (8 entries);
    // duplicate (row, col) pairs are summed only when collapsed to CSR.
    let sum = coo_add(&a, &a);
    assert_eq!(sum.to_csr().values, two_a.to_csr().values);
    println!("coo_add(A, A) collapsed = {:?}", sum.to_csr().values);

    // coo_matvec: A·x with x = [1, 2]  ->  [2·1+1·2, 1·1+3·2] = [4, 7].
    let x = vec![1.0, 2.0];
    let y = coo_matvec(&a, &x);
    println!("coo_matvec(A, [1, 2])  = [{:.1}, {:.1}]", y[0], y[1]);
    assert_eq!(y, vec![4.0, 7.0]);

    // An empty matrix contributes the zero vector of length x.len().
    let empty = Coo::new();
    assert_eq!(coo_matvec(&empty, &x), vec![0.0, 0.0]);
    println!("coo_matvec(empty, x)   = [{:.1}, {:.1}]", 0.0, 0.0);

    println!("\nverified: scale/add/matvec match exact dense arithmetic");
}
