use std::fmt::{Debug, Display};

use clap::ValueEnum;
use num_traits::{Float, Zero};
use numeric_literals::replace_numeric_literals;
use rand_distr::{uniform::SampleUniform, Distribution, Exp, Gamma, Uniform};
use trait_set::trait_set;

trait_set! {
    pub trait BranchLength = Debug + Display + Float + Zero + SampleUniform;
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Distr {
    Uniform,
    Exponential,
    Gamma,
}

pub enum Sampler<T>
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
    pub fn new(v: Distr) -> Self {
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

// #[cfg(test)]
// mod tests {

//     use super::*;

//     #[test]
//     fn sample() {
//         let mut sampler = Sampler::new(Distr::Uniform);
//         let mut rng = rand::thread_rng();
//         for _ in 0..=10 {
//             print!(" {}", sampler.sample(&mut rng))
//         }
//         println!();
//         sampler = Sampler::new(Distr::Exponential);
//         for _ in 0..=10 {
//             print!(" {}", sampler.sample(&mut rng))
//         }
//         println!();
//         sampler = Sampler::new(Distr::Gamma);
//         for _ in 0..=10 {
//             print!(" {}", sampler.sample(&mut rng))
//         }
//         println!();
//         panic!()
//     }
// }
