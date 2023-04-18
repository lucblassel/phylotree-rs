use clap::ValueEnum;
use rand_distr::{Distribution, Exp, Uniform, Gamma};

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Distr {
    Uniform,
    Exponential,
    Gamma,
}

pub enum Sampler {
    Uniform(Uniform<f32>),
    Exponential(Exp<f32>),
    Gamma(Gamma<f32>),
}

impl Sampler {
    pub fn new(v:Distr) -> Self {
        match v {
            Distr::Uniform => Self::Uniform(Uniform::new(0.002, 1.0)),
            Distr::Exponential => Self::Exponential(Exp::new(0.15).unwrap()),
            Distr::Gamma => Self::Gamma(Gamma::new(4.0,1.0).unwrap()),
        }
    }
}

impl Distribution<f32> for Sampler {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> f32 {
        match self {
            Sampler::Uniform(u) => u.sample(rng),
            Sampler::Exponential(e) => e.sample(rng),
            Sampler::Gamma(p) => p.sample(rng),
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn sample() {
        let mut sampler = Sampler::new(Distr::Uniform);
        let mut rng = rand::thread_rng();
        for _ in 0..=10 {
            print!(" {}", sampler.sample(&mut rng))
        }
        println!();
        sampler = Sampler::new(Distr::Exponential);
        for _ in 0..=10 {
            print!(" {}", sampler.sample(&mut rng))
        }
        println!();
        sampler = Sampler::new(Distr::Gamma);
        for _ in 0..=10 {
            print!(" {}", sampler.sample(&mut rng))
        }
        println!();
        panic!()
    }
}