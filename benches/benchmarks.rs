use criterion::BenchmarkId;
use criterion::Criterion;
use criterion::{criterion_group, criterion_main};

use phylotree::{generate_tree, Tree};

fn distance_naive(tree: &Tree) {
    let _matrix = tree.distance_matrix().unwrap();
}

fn distance_recurs(tree: &Tree) {
    let _matrix = tree.distance_matrix_recursive().unwrap();
}

fn from_elem(c: &mut Criterion) {
    let tree: Tree = generate_tree(100, true, phylotree::distr::Distr::Uniform);

    c.bench_with_input(
        BenchmarkId::new("input_uncached", tree.size().unwrap()),
        &tree,
        |b, s| {
            b.iter(|| distance_naive(s));
        },
    );

    c.bench_with_input(
        BenchmarkId::new("input_recursive", tree.size().unwrap()),
        &tree,
        |b, s| {
            b.iter(|| distance_recurs(s));
        },
    );
}

criterion_group!(benches, from_elem);
criterion_main!(benches);
