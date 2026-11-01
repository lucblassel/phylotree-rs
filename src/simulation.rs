//! Seedable topology and branch-length priors.
//!
//! Uniform labeled topologies (PDA), ETE3's deque process and Yule topologies
//! are different distributions. Birth-death uses DendroPy's *first hitting n*
//! convention, restarting extinct realizations and pruning extinct lineages;
//! it does not implement GSA or conditioning on a fixed crown age.
use crate::tree::{Node, NodeId, Tree, TreeError};
use rand::{seq::SliceRandom, Rng};
use rand_distr::{Distribution, Exp, Gamma};
use std::collections::VecDeque;
use thiserror::Error;

/// Failure to validate or realize a simulation prior.
#[derive(Debug, Error)]
pub enum SimulationError {
    /// Invalid prior parameter.
    #[error("Invalid simulation parameter: {0}")]
    InvalidParameter(&'static str),
    /// A resource limit was reached; no alternate prior is silently substituted.
    #[error("Birth-death simulation exceeded its {0} limit")]
    ResourceLimit(&'static str),
    /// Underlying tree operation failed.
    #[error(transparent)]
    Tree(#[from] TreeError),
}

/// Independent edge-length prior, or the birth-death process's elapsed times.
#[derive(Debug, Clone, Copy)]
pub enum BranchLengths {
    /// Exponential distribution parameterized by its mean, not its rate.
    Exponential {
        /// Mean edge length, strictly positive.
        mean: f64,
    },
    /// Gamma distribution parameterized by shape and scale (mean = shape * scale).
    Gamma {
        /// Positive shape.
        shape: f64,
        /// Positive scale, not rate.
        scale: f64,
    },
    /// A fixed nonnegative edge length (also useful for topology-only studies).
    Fixed(f64),
    /// Keep elapsed birth-death times, including the two zero-length newborn edges.
    ProcessTime,
}
impl Default for BranchLengths {
    fn default() -> Self {
        Self::Exponential { mean: 0.1 }
    }
}

/// Constant-rate birth-death process and explicit resource limits.
#[derive(Debug, Clone, Copy)]
pub struct BirthDeath {
    /// Per-lineage birth rate (strictly positive).
    pub birth_rate: f64,
    /// Per-lineage death rate (nonnegative).
    pub death_rate: f64,
    /// Maximum total events, including extinct restarts. Reaching it returns an error.
    pub max_events: usize,
    /// Maximum extinct restarts before returning an error.
    pub max_restarts: usize,
}
impl Default for BirthDeath {
    fn default() -> Self {
        Self {
            birth_rate: 1.0,
            death_rate: 0.5,
            max_events: 10_000_000,
            max_restarts: 100_000,
        }
    }
}

/// The distribution on topology; the root is retained for the rooted processes.
#[derive(Debug, Clone, Copy)]
pub enum Topology {
    /// Uniform on (2n-5)!! unrooted labeled binary topologies.
    UniformUnrooted,
    /// Uniform on (2n-3)!! rooted labeled binary topologies, including stem insertions.
    UniformRooted,
    /// ETE3 3.1.3 Tree.populate: pop either end of a leaf deque with equal probability.
    Ete3,
    /// Yule/Harding topology: choose a uniformly random extant leaf to split.
    Yule,
    /// Start with one lineage; stop immediately on first reaching n extant tips.
    BirthDeath(BirthDeath),
}

/// Generate a tree with caller-owned RNG state and unique labels T1..Tn.
///
/// `ProcessTime` is valid only for `BirthDeath`. Other branch priors replace
/// all non-root edges *after* topology generation. The root stem is omitted.
///
/// ```
/// use phylotree::simulation::{simulate, Topology, BranchLengths};
/// use rand::{SeedableRng, rngs::StdRng};
/// let mut rng = StdRng::seed_from_u64(42);
/// let t = simulate(20, Topology::UniformUnrooted,
///                  BranchLengths::Gamma { shape: 2.0, scale: 0.05 }, &mut rng).unwrap();
/// assert_eq!(t.n_leaves(), 20);
/// ```
pub fn simulate<R: Rng + ?Sized>(
    tips: usize,
    topology: Topology,
    lengths: BranchLengths,
    rng: &mut R,
) -> Result<Tree, SimulationError> {
    if tips < 2
        || tips > usize::MAX / 4
        || (tips < 3 && matches!(topology, Topology::UniformUnrooted))
    {
        return Err(SimulationError::InvalidParameter("tip count"));
    }
    if matches!(lengths, BranchLengths::ProcessTime) && !matches!(topology, Topology::BirthDeath(_))
    {
        return Err(SimulationError::InvalidParameter(
            "process times require a birth-death process",
        ));
    }
    let sampler = LengthSampler::new(lengths)?;
    let mut tree = match topology {
        Topology::UniformUnrooted => uniform_unrooted(tips, rng),
        Topology::UniformRooted => uniform_rooted(tips, rng),
        Topology::Ete3 | Topology::Yule => {
            leaf_splitting(tips, matches!(topology, Topology::Ete3), rng)
        }
        Topology::BirthDeath(prior) => birth_death(tips, prior, rng)?,
    };
    if !matches!(lengths, BranchLengths::ProcessTime) {
        for id in preorder(&tree)? {
            let parent = tree.get(id)?.parent;
            if let Some(parent) = parent {
                let length = sampler.sample(rng);
                if !length.is_finite() {
                    return Err(SimulationError::InvalidParameter(
                        "sampled nonfinite branch length",
                    ));
                }
                tree.get_mut(id)?.parent_edge = Some(length);
                tree.get_mut(parent)?.set_child_edge(id, Some(length));
            }
        }
    }
    Ok(tree)
}

enum LengthSampler {
    Exp(Exp<f64>),
    Gamma(Gamma<f64>),
    Fixed(f64),
    Process,
}
impl LengthSampler {
    fn new(prior: BranchLengths) -> Result<Self, SimulationError> {
        let positive = |x: f64| x.is_finite() && x > 0.0;
        Ok(match prior {
            BranchLengths::Exponential { mean } if positive(mean) && positive(1.0 / mean) => {
                Self::Exp(Exp::new(1.0 / mean).unwrap())
            }
            BranchLengths::Gamma { shape, scale } if positive(shape) && positive(scale) => {
                Self::Gamma(Gamma::new(shape, scale).unwrap())
            }
            BranchLengths::Fixed(x) if x.is_finite() && x >= 0.0 => Self::Fixed(x),
            BranchLengths::ProcessTime => Self::Process,
            _ => return Err(SimulationError::InvalidParameter("branch-length prior")),
        })
    }
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        match self {
            Self::Exp(d) => d.sample(rng),
            Self::Gamma(d) => d.sample(rng),
            Self::Fixed(x) => *x,
            Self::Process => unreachable!(),
        }
    }
}

// Construct both parent and child edge representations once, without repeated
// insertion/cache work. Only connected live nodes are supplied by the samplers.
fn build(
    children: Vec<Vec<usize>>,
    lengths: Vec<f64>,
    root: usize,
    labels: Vec<Option<String>>,
) -> Tree {
    let mut nodes: Vec<_> = labels
        .into_iter()
        .enumerate()
        .map(|(i, name)| {
            let mut n = Node::new();
            n.id = i.into();
            n.name = name;
            n
        })
        .collect();
    let mut stack = vec![root];
    while let Some(p) = stack.pop() {
        for &c in &children[p] {
            nodes[c].set_parent(p.into(), Some(lengths[c]));
            let depth = nodes[p].get_depth() + 1;
            nodes[c].set_depth(depth);
            nodes[p].add_child(c.into(), Some(lengths[c]));
            stack.push(c);
        }
    }
    Tree::from_nodes(nodes)
}
fn random_labels<R: Rng + ?Sized>(children: &[Vec<usize>], rng: &mut R) -> Vec<Option<String>> {
    let mut names: Vec<_> = (1..=children.iter().filter(|c| c.is_empty()).count())
        .map(|i| format!("T{i}"))
        .collect();
    names.shuffle(rng);
    children
        .iter()
        .map(|c| if c.is_empty() { names.pop() } else { None })
        .collect()
}
fn uniform_unrooted<R: Rng + ?Sized>(tips: usize, rng: &mut R) -> Tree {
    let mut edges = vec![(0, tips), (1, tips), (2, tips)];
    for tip in 3..tips {
        let i = rng.gen_range(0..edges.len());
        let (a, b) = edges[i];
        let internal = tips + tip - 2;
        edges[i] = (a, internal);
        edges.push((internal, b));
        edges.push((internal, tip));
    }
    let count = 2 * tips - 2;
    let mut adj = vec![Vec::new(); count];
    for (a, b) in edges {
        adj[a].push(b);
        adj[b].push(a);
    }
    let mut children = vec![Vec::new(); count];
    let mut stack = vec![(tips, usize::MAX)];
    while let Some((p, dad)) = stack.pop() {
        for &c in &adj[p] {
            if c != dad {
                children[p].push(c);
                stack.push((c, p));
            }
        }
    }
    let names = (0..count)
        .map(|i| (i < tips).then(|| format!("T{}", i + 1)))
        .collect();
    build(children, vec![1.0; count], tips, names)
}
fn uniform_rooted<R: Rng + ?Sized>(tips: usize, rng: &mut R) -> Tree {
    let mut children = vec![vec![1, 2], vec![], vec![]];
    let mut parents = vec![None, Some(0), Some(0)];
    let mut root = 0;
    for _ in 2..tips {
        // Every node identifies an edge above it; root identifies the stem edge.
        let old = rng.gen_range(0..children.len());
        let internal = children.len();
        let tip = internal + 1;
        let parent = parents[old];
        if let Some(p) = parent {
            let pos = children[p].iter().position(|&i| i == old).unwrap();
            children[p][pos] = internal;
        } else {
            root = internal;
        }
        children.push(vec![old, tip]);
        children.push(vec![]);
        parents.push(parent);
        parents.push(Some(internal));
        parents[old] = Some(internal);
    }
    let names = random_labels(&children, rng);
    let count = children.len();
    build(children, vec![1.0; count], root, names)
}

fn leaf_splitting<R: Rng + ?Sized>(tips: usize, ete: bool, rng: &mut R) -> Tree {
    let mut children = vec![vec![]];
    let mut queue = VecDeque::from([0]);
    let mut active = vec![0];
    for _ in 1..tips {
        let p = if ete {
            if rng.gen_bool(0.5) {
                queue.pop_front().unwrap()
            } else {
                queue.pop_back().unwrap()
            }
        } else {
            let i = rng.gen_range(0..active.len());
            active.swap_remove(i)
        };
        let a = children.len();
        children.push(vec![]);
        children.push(vec![]);
        children[p] = vec![a, a + 1];
        if ete {
            queue.extend([a, a + 1]);
        } else {
            active.extend([a, a + 1]);
        }
    }
    let names = random_labels(&children, rng);
    let count = children.len();
    build(children, vec![1.0; count], 0, names)
}

fn birth_death<R: Rng + ?Sized>(
    tips: usize,
    prior: BirthDeath,
    rng: &mut R,
) -> Result<Tree, SimulationError> {
    let total = prior.birth_rate + prior.death_rate;
    if !prior.birth_rate.is_finite()
        || prior.birth_rate <= 0.0
        || !prior.death_rate.is_finite()
        || prior.death_rate < 0.0
        || !total.is_finite()
        || prior.max_events == 0
        || prior.max_restarts == 0
    {
        return Err(SimulationError::InvalidParameter(
            "birth/death rates or limits",
        ));
    }
    #[derive(Default)]
    struct Lineage {
        children: Option<[usize; 2]>,
        end: Option<f64>,
        alive: bool,
    }
    let mut events = 0;
    for _ in 0..prior.max_restarts {
        let mut nodes = vec![Lineage {
            alive: true,
            ..Default::default()
        }];
        let mut active = vec![0];
        let mut time = 0.0;
        while !active.is_empty() && active.len() < tips {
            events += 1;
            if events > prior.max_events {
                return Err(SimulationError::ResourceLimit("event"));
            }
            let rate = total * active.len() as f64;
            if !rate.is_finite() {
                return Err(SimulationError::InvalidParameter("total event rate"));
            }
            time += Exp::new(rate)
                .map_err(|_| SimulationError::InvalidParameter("event rate"))?
                .sample(rng);
            if !time.is_finite() {
                return Err(SimulationError::InvalidParameter("nonfinite event time"));
            }
            let i = rng.gen_range(0..active.len());
            let node = active.swap_remove(i);
            nodes[node].alive = false;
            nodes[node].end = Some(time);
            if rng.gen_bool(prior.birth_rate / total) {
                let a = nodes.len();
                nodes[node].children = Some([a, a + 1]);
                nodes.push(Lineage {
                    alive: true,
                    ..Default::default()
                });
                nodes.push(Lineage {
                    alive: true,
                    ..Default::default()
                });
                active.extend([a, a + 1]);
            }
        }
        if active.is_empty() {
            continue;
        }
        let mut survivors = vec![0; nodes.len()];
        for i in (0..nodes.len()).rev() {
            survivors[i] = if let Some([a, b]) = nodes[i].children {
                survivors[a] + survivors[b]
            } else {
                usize::from(nodes[i].alive)
            };
        }
        let mut children: Vec<Vec<usize>> = Vec::new();
        let mut lengths = Vec::new();
        let mut stack: Vec<(usize, Option<usize>, f64)> = vec![(0, None, 0.0)];
        while let Some((old, parent, parent_time)) = stack.pop() {
            if survivors[old] == 0 {
                continue;
            }
            let live: Vec<_> = nodes[old]
                .children
                .into_iter()
                .flatten()
                .filter(|&c| survivors[c] > 0)
                .collect();
            if live.len() == 1 {
                stack.push((live[0], parent, parent_time));
                continue;
            }
            let id = children.len();
            children.push(vec![]);
            let end = nodes[old].end.unwrap_or(time);
            lengths.push(if parent.is_some() {
                end - parent_time
            } else {
                0.0
            });
            if let Some(p) = parent {
                children[p].push(id);
            }
            for &c in live.iter().rev() {
                stack.push((c, Some(id), end));
            }
        }
        let names = random_labels(&children, rng);
        return Ok(build(children, lengths, 0, names));
    }
    Err(SimulationError::ResourceLimit("extinct-restart"))
}

fn preorder(tree: &Tree) -> Result<Vec<NodeId>, SimulationError> {
    let mut pending = vec![tree.root()?];
    let mut order = Vec::new();
    while let Some(id) = pending.pop() {
        order.push(id);
        pending.extend(tree.get(id)?.children.iter().rev().copied());
    }
    Ok(order)
}
