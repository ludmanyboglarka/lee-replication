//! TRANSLATION OF RDMC_SOA.m
//! = RDMC trial simulator, first passage rule added: ? 
//!
//! @inputs: A, k, v_a, v_c, c
//! @output: an array of 1500 positions: position i is the state i at i ms. 
//! function knows nothing about boundaries -> trial ending and correct/incorrect not implemented
//!
//! w_a(t) = A * exp(-k * t)
//! v(t)   = w_a(t) * v_a + (1 - w_a(t)) * v_c
//! X(t+dt) = X(t) + v(t) * dt + SIGMA * sqrt(dt) * N(0, 1)
//!
//! Stops at the first boundary crossed (+a correct, -a error). No path is stored.
//! Non-decision time is added by the caller.

use rand::Rng;
use rand_distr::{Distribution, StandardNormal};

// fixed parameters: 1500 steps
// pub: other files in this project can use it
pub const SIGMA: f64 = 0.1; 
pub const DT: f64 = 0.001; // f64: floating-point num, decimal values
pub const N_STEPS: u32 = 1500; // u32: non-negative whole number

#[derive(Debug, Clone, Copy)] 
pub struct TrialParams {
    pub a_weight: f64, // A
    pub k: f64,
    pub v_a: f64, // +da congruent, -da incongruent
    pub v_c: f64,
    pub start: f64, // c
    pub boundary: f64, // a
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    Correct { step: u32 }, // step = decision time in ms
    Error { step: u32 },
    Timeout,
}

// fn simulates one trial: particle drifts and jitters until it hits a boundary or time runs out
pub fn simulate_trial<R: Rng + ?Sized>(p: &TrialParams, rng: &mut R) -> Outcome {
    let noise_sd = SIGMA * DT.sqrt();
    let mut x = p.start;

    for step in 1..=N_STEPS {
        let t = f64::from(step) * DT;
        // drift
        let w_a = p.a_weight * (-p.k * t).exp();
        let drift = w_a * p.v_a + (1.0 - w_a) * p.v_c;

        let z: f64 = StandardNormal.sample(rng);
        x += drift * DT + noise_sd * z; // The evidence moves by a drift part (direction and speed) plus a random part

        // Checking the boundaries
        if x >= p.boundary {
            return Outcome::Correct { step };
        }
        if x <= -p.boundary {
            return Outcome::Error { step };
        }
    }
    Outcome::Timeout
}

// test code: only runs with cargo test -> job is to test whether the simulation is correct
#[cfg(test)]
mod tests {
    use super::*; // "bring in everything from the file above me"
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    // A = 0, c = 0: compare with closed-form accuracy and mean decision time.
    // a is shifted by 0.5826*sigma*sqrt(dt) to account for the 1 ms grid.
    // so: if the autom. process is switched off (A = 0), and the start is in the middle (c =
    // 0), the model is a plain ddm -> the formulas give the exact acc. and mean decision time.
    //   many trials are simulated, compute the same 2 numbers from the sim., and check if they
    //   agree with the formulas -> if so, the sim. is almost certainly right. 
    #[test]
    fn matches_closed_form() {
        for &(v, a) in &[(0.3, 0.1), (0.4, 0.1)] {
            let p = TrialParams { a_weight: 0.0, k: 0.0, v_a: 0.0, v_c: v, start: 0.0, boundary: a };
            let mut rng = StdRng::seed_from_u64(1);
            let (mut n_c, mut n_e, mut sum_c) = (0u32, 0u32, 0.0f64);
            for _ in 0..500_000 {
                match simulate_trial(&p, &mut rng) {
                    Outcome::Correct { step } => { n_c += 1; sum_c += f64::from(step); }
                    Outcome::Error { .. } => n_e += 1,
                    Outcome::Timeout => {}
                }
            }
            let acc = f64::from(n_c) / f64::from(n_c + n_e);
            let mean_c = sum_c / f64::from(n_c);

            let s2 = SIGMA * SIGMA;
            let a_eff = a + 0.5826 * SIGMA * DT.sqrt();
            let acc_exp = 1.0 / (1.0 + (-2.0 * v * a_eff / s2).exp());
            let t_exp = (a_eff / v) * (v * a_eff / s2).tanh() * 1000.0;

            println!("v={v}: acc {acc:.4} (exp {acc_exp:.4}), mean RT {mean_c:.1} (exp {t_exp:.1})");
            assert!((acc - acc_exp).abs() < 0.003);
            assert!((mean_c - t_exp).abs() / t_exp < 0.01);
        }
    }

    #[test]
    fn same_seed_same_result() {
        let p = TrialParams { a_weight: 0.5, k: 8.0, v_a: 0.3, v_c: 0.2, start: 0.0, boundary: 0.1 };
        let run = || {
            let mut rng = StdRng::seed_from_u64(42);
            (0..100).map(|_| simulate_trial(&p, &mut rng)).collect::<Vec<_>>()
        };
        assert_eq!(run(), run());
    }
}
