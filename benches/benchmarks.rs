use criterion::{criterion_group, criterion_main, PlotConfiguration};
use criterion::{BenchmarkId, Criterion};

use itertools::Itertools;
use ndarray::{Array2, Axis};
use ndarray_rand::rand_distr::Uniform as UniformND;
use ndarray_rand::RandomExt;
use phylotree::distance::DistanceMatrix;
use phylotree::{distr::Distr::Uniform, generate_tree};

/// Measure how distance matrix extraction scales with tree size
fn dm_vs_treesize(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(criterion::AxisScale::Logarithmic);

    let mut group = c.benchmark_group("dm_vs_treesize");
    group.plot_config(plot_config);

    for size in [10, 20, 40, 100, 500, 1000, 2000, 5000, 10000].iter() {
        group.bench_with_input(BenchmarkId::new("PhyloDM", size), size, |bencher, size| {
            bencher.iter(|| {
                let tree = generate_tree(*size, true, Uniform).unwrap();
                tree.distance_matrix()
            })
        });
        if *size <= 1000 {
            group.bench_with_input(
                BenchmarkId::new("Recursive", size),
                size,
                |bencher, size| {
                    bencher.iter(|| {
                        let tree = generate_tree(*size, true, Uniform).unwrap();
                        tree.distance_matrix_recursive()
                    })
                },
            );
        }
    }

    group.finish();
}

/// Benchmark newick parsing
fn newick_parsing(c: &mut Criterion) {
    let mut group = c.benchmark_group("newick_parsing");
    for size in [10, 20, 40, 100, 500, 1000, 2000, 5000, 10000, 20000].iter() {
        let newick = generate_tree(*size, true, Uniform)
            .unwrap()
            .to_newick()
            .unwrap();
        group.bench_with_input(BenchmarkId::from_parameter(*size), size, |bencher, _| {
            bencher.iter(|| {
                let _ = phylotree::tree::Tree::from_newick(&newick);
            })
        });
    }
}

fn generate_phylip(size: usize) -> String {
    // Generate random matrix
    let base = Array2::random((size, size), UniformND::new(0., 1.));
    let base = base.view();
    // Make it symmetric
    let mut sym = &base + &base.t();
    sym.diag_mut().fill(0.);

    let s = sym
        .axis_iter(Axis(0))
        .zip(0..size)
        .map(|(row, i)| {
            let vals = row.iter().map(|v| v.to_string()).join("  ");
            format!("T{i}    {vals}")
        })
        .join("\n");

    format!("{size}\n{s}\n")
}

fn phylip_parsing(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(criterion::AxisScale::Logarithmic);
    let mut group = c.benchmark_group("phylip_parsing");
    group.plot_config(plot_config);

    for size in [10, 20, 40, 100, 500, 1000, 2000, 5000, 10000].iter() {
        let phylip = generate_phylip(*size);
        group.bench_with_input(BenchmarkId::new("Lower", size), size, |bencher, _| {
            bencher.iter(|| {
                let _ = DistanceMatrix::<f64>::from_phylip_tril(&phylip);
            })
        });
        if *size <= 1000 {
            group.bench_with_input(BenchmarkId::new("Strict", size), size, |bencher, _| {
                bencher.iter(|| {
                    let _ = DistanceMatrix::<f64>::from_phylip_strict(&phylip, true);
                })
            });
        }
    }
}

criterion_group!(benches, dm_vs_treesize, newick_parsing, phylip_parsing);
criterion_main!(benches);
