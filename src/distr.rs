//! Distributions to generate branch lengths inb random trees
//!

use std::fmt::{Debug, Display};

use clap::ValueEnum;
use num_traits::{Float, Zero};
use numeric_literals::replace_numeric_literals;
use rand_distr::{uniform::SampleUniform, Distribution, Exp, Gamma, Uniform};
use trait_set::trait_set;

trait_set! {
    /// Trait describing objects that can be used as branch lengths
    /// in phylogenetic trees.
    pub trait BranchLength = Debug + Display + Float + Zero + SampleUniform;
}

/// Available branch length distributions
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Distr {
    /// A [uniform](https://en.wikipedia.org/wiki/Continuous_uniform_distribution)
    /// distribution over $[0.002, 1.0)$
    Uniform,
    /// An [exponential](https://en.wikipedia.org/wiki/Exponential_distribution)
    /// distribution with rate $\lambda=0.15$
    Exponential,
    /// A [gamma](https://en.wikipedia.org/wiki/Gamma_distribution) distribution
    /// with a shape $k=4$ and scale $\theta=1.0$.
    Gamma,
}

pub(crate) enum Sampler<T>
where
    T: BranchLength,
    rand_distr::StandardNormal: rand_distr::Distribution<T>,
    rand_distr::Exp1: rand_distr::Distribution<T>,
    rand_distr::Open01: rand_distr::Distribution<T>,
{
    Uniform(Uniform<T>),
    Exponential(Exp<T>),
    Gamma(Gamma<T>),
}

impl<T> Sampler<T>
where
    T: BranchLength,
    rand_distr::StandardNormal: rand_distr::Distribution<T>,
    rand_distr::Exp1: rand_distr::Distribution<T>,
    rand_distr::Open01: rand_distr::Distribution<T>,
{
    #[replace_numeric_literals(T::from(literal).unwrap())]
    pub(crate) fn new(v: Distr) -> Self {
        match v {
            Distr::Uniform => Self::Uniform(Uniform::<T>::new(0.002, 1.0)),
            Distr::Exponential => Self::Exponential(Exp::new(0.15).unwrap()),
            Distr::Gamma => Self::Gamma(Gamma::new(4.0, 1.0).unwrap()),
        }
    }
}

impl<T> Distribution<T> for Sampler<T>
where
    T: BranchLength,
    rand_distr::StandardNormal: rand_distr::Distribution<T>,
    rand_distr::Exp1: rand_distr::Distribution<T>,
    rand_distr::Open01: rand_distr::Distribution<T>,
{
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> T {
        match self {
            Sampler::Uniform(u) => u.sample(rng),
            Sampler::Exponential(e) => e.sample(rng),
            Sampler::Gamma(p) => p.sample(rng),
        }
    }
}
