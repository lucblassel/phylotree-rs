use phylotree::{
    simulation::{simulate, BirthDeath, BranchLengths, Topology},
    tree::Tree,
};
use rand::{rngs::StdRng, SeedableRng};
use std::collections::{HashMap, HashSet};
fn signature(tree: &Tree, rooted: bool) -> String {
    let root = tree.root().unwrap();
    let mut sets: HashMap<phylotree::tree::NodeId, Vec<String>> = HashMap::new();
    let mut splits = Vec::new();
    let mut all: Vec<_> = tree.leaf_names().into_iter().flatten().collect();
    all.sort();
    for id in tree.postorder(root).unwrap() {
        let node = tree.get(id).unwrap();
        let mut names: Vec<String> = if node.children().is_empty() {
            vec![node.name().unwrap().to_owned()]
        } else {
            node.children()
                .iter()
                .flat_map(|c| sets[c].clone())
                .collect()
        };
        names.sort();
        if id != root && names.len() > 1 && names.len() < all.len() - (usize::from(!rooted)) {
            let complement: Vec<_> = all.iter().filter(|n| !names.contains(n)).cloned().collect();
            let canonical = if !rooted
                && (complement.len() < names.len()
                    || (complement.len() == names.len() && complement < names))
            {
                complement
            } else {
                names.clone()
            };
            splits.push(canonical.join(","));
        }
        sets.insert(id, names);
    }
    splits.sort();
    splits.join("|")
}
#[test]
fn uniform_labeled_topology_frequencies_and_distinct_ete_yule_priors() {
    let mut rng = StdRng::seed_from_u64(101);
    for (top, tips, rooted) in [
        (Topology::UniformUnrooted, 5, false),
        (Topology::UniformRooted, 4, true),
    ] {
        let mut counts = HashMap::new();
        for _ in 0..30_000 {
            let t = simulate(tips, top, BranchLengths::Fixed(1.0), &mut rng).unwrap();
            *counts.entry(signature(&t, rooted)).or_insert(0_usize) += 1;
        }
        assert_eq!(counts.len(), 15);
        for c in counts.values() {
            assert!((*c as f64 - 2000.0).abs() < 220.0, "{counts:?}");
        }
    }
    for (top, expected) in [
        (Topology::UniformRooted, 0.2),
        (Topology::Ete3, 0.5),
        (Topology::Yule, 1.0 / 3.0),
    ] {
        let mut balanced = 0;
        for _ in 0..20_000 {
            let t = simulate(4, top, BranchLengths::Fixed(1.0), &mut rng).unwrap();
            let r = t.get(t.root().unwrap()).unwrap();
            if r.children()
                .iter()
                .all(|&c| t.get(c).unwrap().children().len() == 2)
            {
                balanced += 1;
            }
        }
        let observed = balanced as f64 / 20_000.0;
        assert!((observed - expected).abs() < 0.015, "{top:?}: {observed}");
    }
}
#[test]
fn branch_moments_reproducibility_and_diameter() {
    for (prior, mean, variance) in [
        (BranchLengths::Exponential { mean: 0.1 }, 0.1, 0.01),
        (
            BranchLengths::Gamma {
                shape: 2.0,
                scale: 0.05,
            },
            0.1,
            0.005,
        ),
    ] {
        let mut rng = StdRng::seed_from_u64(77);
        let mut values = Vec::new();
        for _ in 0..2500 {
            let mut t = simulate(10, Topology::UniformUnrooted, prior, &mut rng).unwrap();
            for id in t.preorder(t.root().unwrap()).unwrap() {
                let n = t.get(id).unwrap();
                if n.parent().is_some() {
                    values.push(n.parent_edge().unwrap());
                }
            }
            if values.len() < 200 {
                assert!((pairwise_diameter(&t) - t.diameter().unwrap()).abs() < 1e-12);
                t.rescale_diameter(3.0).unwrap();
                assert!((t.diameter().unwrap() - 3.0).abs() < 1e-12);
            }
        }
        let m = values.iter().sum::<f64>() / values.len() as f64;
        let v = values.iter().map(|x| (x - m).powi(2)).sum::<f64>() / values.len() as f64;
        assert!((m - mean).abs() < 0.003);
        assert!((v - variance).abs() < 0.001);
    }
    for top in [
        Topology::UniformUnrooted,
        Topology::UniformRooted,
        Topology::Ete3,
        Topology::Yule,
        Topology::BirthDeath(BirthDeath::default()),
    ] {
        let mut a = StdRng::seed_from_u64(42);
        let mut b = StdRng::seed_from_u64(42);
        assert_eq!(
            simulate(20, top, BranchLengths::default(), &mut a)
                .unwrap()
                .to_newick()
                .unwrap(),
            simulate(20, top, BranchLengths::default(), &mut b)
                .unwrap()
                .to_newick()
                .unwrap()
        );
    }
}
#[test]
fn birth_death_first_hitting_time_and_pruning() {
    let mut rng = StdRng::seed_from_u64(981);
    let mut heights = Vec::new();
    for death in [0.0, 0.5] {
        for _ in 0..3000 {
            let t = simulate(
                10,
                Topology::BirthDeath(BirthDeath {
                    death_rate: death,
                    ..Default::default()
                }),
                BranchLengths::ProcessTime,
                &mut rng,
            )
            .unwrap();
            let root = t.root().unwrap();
            let leaves = t.leaves();
            assert_eq!(leaves.len(), 10);
            assert_eq!(t.preorder(root).unwrap().len(), 19);
            assert_eq!(t.get(root).unwrap().children().len(), 2);
            assert_eq!(
                t.leaf_names()
                    .into_iter()
                    .flatten()
                    .collect::<HashSet<_>>()
                    .len(),
                10
            );
            let h: Vec<_> = leaves
                .iter()
                .map(|&id| t.get_distance(root, id).unwrap().0.unwrap())
                .collect();
            assert!(h.iter().all(|x| (x - h[0]).abs() < 1e-12));
            assert_eq!(
                leaves
                    .iter()
                    .filter(|&&id| t.get(id).unwrap().parent_edge().unwrap() == 0.0)
                    .count(),
                2
            );
            if death == 0.0 {
                heights.push(h[0]);
            }
        }
    }
    let expectation = (2..10).map(|k| 1.0 / k as f64).sum::<f64>();
    let observed = heights.iter().sum::<f64>() / heights.len() as f64;
    assert!(
        (observed - expectation).abs() < 0.05,
        "{observed} vs {expectation}"
    );
}
#[test]
fn invalid_parameters_and_limits_do_not_change_the_prior() {
    let mut rng = StdRng::seed_from_u64(1);
    assert!(simulate(1, Topology::Ete3, BranchLengths::default(), &mut rng).is_err());
    for prior in [
        BranchLengths::Exponential { mean: 0.0 },
        BranchLengths::Gamma {
            shape: -1.0,
            scale: 1.0,
        },
        BranchLengths::ProcessTime,
    ] {
        assert!(simulate(10, Topology::Ete3, prior, &mut rng).is_err());
    }
    assert!(simulate(
        10,
        Topology::BirthDeath(BirthDeath {
            max_events: 1,
            ..Default::default()
        }),
        BranchLengths::ProcessTime,
        &mut rng
    )
    .is_err());
    let mut zero = simulate(10, Topology::Ete3, BranchLengths::Fixed(0.0), &mut rng).unwrap();
    assert!(zero.rescale_diameter(1.0).is_err());
}

// Independent oracle retains the previous pairwise diameter definition.
fn pairwise_diameter(tree: &Tree) -> f64 {
    let leaves = tree.leaves();
    let mut distances = Vec::new();
    for (i, &a) in leaves.iter().enumerate() {
        for &b in &leaves[i + 1..] {
            let (weighted, hops) = tree.get_distance(a, b).unwrap();
            distances.push(weighted.unwrap_or(hops as f64));
        }
    }
    distances.into_iter().reduce(f64::max).unwrap()
}

#[test]
fn linear_diameter_preserves_pairwise_semantics() {
    use rand::Rng;
    let mut rng = StdRng::seed_from_u64(128);
    for tips in [2, 3, 10, 50] {
        for _ in 0..30 {
            let original =
                simulate(tips, Topology::Yule, BranchLengths::default(), &mut rng).unwrap();
            for mode in 0..4 {
                let mut tree = original.clone();
                for id in tree.preorder(tree.root().unwrap()).unwrap() {
                    let node = tree.get_mut(id).unwrap();
                    let edge = match mode {
                        1 => None,
                        2 if rng.gen_bool(0.4) => None,
                        3 => Some(rng.gen_range(-3.0..3.0)),
                        _ => node.parent_edge(),
                    };

                    if let Some(parent) = node.parent() {
                        tree.set_edge(parent, id, edge).unwrap()
                    }
                }
                assert!((tree.diameter().unwrap() - pairwise_diameter(&tree)).abs() < 1e-12);
            }
        }
    }
    for newick in [
        "(A:100,B,C);",
        "((A:2):1,B:4,C:5):100;",
        "(A:-3,B:-2,C:-1);",
    ] {
        let tree = Tree::from_newick(newick).unwrap();
        assert_eq!(tree.diameter().unwrap(), pairwise_diameter(&tree));
    }
    assert!(Tree::new().diameter().is_err());
    assert!(Tree::from_newick("(A:1);").unwrap().diameter().is_err());
}

#[test]
fn diameter_rescaling_rejects_invalid_inputs_without_mutation() {
    for newick in ["(A:1,B:2);", "(A,B);", "(A:-1,B:2);", "(A:0,B:0);"] {
        for target in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut tree = Tree::from_newick(newick).unwrap();
            let before = tree.to_newick().unwrap();
            assert!(tree.rescale_diameter(target).is_err());
            assert_eq!(before, tree.to_newick().unwrap());
        }
    }
    for newick in ["(A,B);", "(A:-1,B:2);", "(A:0,B:0);"] {
        let mut tree = Tree::from_newick(newick).unwrap();
        let before = tree.to_newick().unwrap();
        assert!(tree.rescale_diameter(3.0).is_err());
        assert_eq!(before, tree.to_newick().unwrap());
    }
    let mut tree = Tree::from_newick("((A:1):2,B:3,C:4):99;").unwrap();
    tree.rescale_diameter(14.0).unwrap();
    assert_eq!(tree.diameter().unwrap(), 14.0);
    assert_eq!(
        tree.get(tree.root().unwrap()).unwrap().parent_edge(),
        Some(198.0)
    );
}

#[test]
fn diameter_handles_deep_trees_without_recursive_traversal() {
    use phylotree::tree::Node;
    let mut tree = Tree::new();
    let root = tree.add_root(Node::new()).unwrap();
    tree.add_child(Node::new(), root, Some(1.0)).unwrap();
    let mut spine = root;
    for _ in 0..10_000 {
        spine = tree.add_child(Node::new(), spine, Some(1.0)).unwrap();
    }
    assert_eq!(tree.diameter().unwrap(), 10_001.0);
    tree.rescale_diameter(20_002.0).unwrap();
    assert_eq!(tree.diameter().unwrap(), 20_002.0);
}
