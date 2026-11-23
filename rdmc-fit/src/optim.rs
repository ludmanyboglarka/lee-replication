//! Bounded Nelder-Mead 
//! maps between "external" parameters (real values, inside their bounds)
//!   and "internal" ones (unbounded, what the simplex moves in), as in `fminsearchbnd`:
//!     two-sided [lo, hi], lower bound only, upper bound only, no bounds
use anyhow::{ensure, Result};

const ZERO_STEP: f64 = 0.00025;

// ---------------------------------------------------------------------------------
// Bounds
// ---------------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct Bounds {
    limits: Vec<(f64, f64)>,
}

fn to_external_1(u: f64, lo: f64, hi: f64) -> f64 {
    match (lo.is_finite(), hi.is_finite()) {
        (true, true) => lo + 0.5 * (u.sin() + 1.0) * (hi - lo),
        (true, false) => lo + u * u,
        (false, true) => hi - u * u,
        (false, false) => u,
    }
}

fn to_internal_1(x: f64, lo: f64, hi: f64) -> f64 {
    match (lo.is_finite(), hi.is_finite()) {
        (true, true) => (2.0 * (x - lo) / (hi - lo) - 1.0).clamp(-1.0, 1.0).asin(),
        (true, false) => (x - lo).sqrt(),
        (false, true) => (hi - x).sqrt(),
        (false, false) => x,
    }
}

impl Bounds {
    /// `limits[i] = (lower, upper)`; use `f64::NEG_INFINITY` / `f64::INFINITY` for "no bound".
    pub fn new(limits: Vec<(f64, f64)>) -> Result<Self> {
        for (i, &(lo, hi)) in limits.iter().enumerate() {
            ensure!(
                !lo.is_nan() && !hi.is_nan() && lo < hi,
                "bounds of parameter {i} must satisfy lower < upper, got ({lo}, {hi})"
            );
        }
        Ok(Self { limits })
    }

    pub fn unbounded(n: usize) -> Self {
        Self { limits: vec![(f64::NEG_INFINITY, f64::INFINITY); n] }
    }

    pub fn len(&self) -> usize {
        self.limits.len()
    }

    pub fn is_empty(&self) -> bool {
        self.limits.is_empty()
    }

    /// Internal (unbounded) -> external (inside the bounds)
    pub fn to_external(&self, u: &[f64]) -> Vec<f64> {
        u.iter()
            .zip(&self.limits)
            .map(|(&u, &(lo, hi))| to_external_1(u, lo, hi))
            .collect()
    }

    /// External -> internal
    pub fn to_internal(&self, x: &[f64]) -> Result<Vec<f64>> {
        ensure!(x.len() == self.limits.len(), "expected {} values, got {}", self.limits.len(), x.len());
        let mut out = Vec::with_capacity(x.len());
        for (i, (&x, &(lo, hi))) in x.iter().zip(&self.limits).enumerate() {
            ensure!(
                x.is_finite() && x >= lo && x <= hi,
                "parameter {i} = {x} is outside its bounds [{lo}, {hi}]"
            );
            out.push(to_internal_1(x, lo, hi));
        }
        Ok(out)
    }
}

// ---------------------------------------------------------------------------------
// Nelder-Mead
// ---------------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopReason {
    Converged,
    MaxIter,
    MaxEvals,
}

#[derive(Debug, Clone)]
pub struct Solution {
    /// Best parameters (external, inside the bounds)
    pub x: Vec<f64>,
    pub f: f64,
    pub evals: usize,
    pub iters: usize,
    pub stop: StopReason,
}

#[derive(Debug, Clone, Copy)]
pub struct NelderMead {
    pub max_iter: usize,
    /// Checked once per iteration, so it can be exceeded by a few evaluations.
    pub max_evals: usize,
    /// Simplex size: largest coordinate distance to the best vertex, relative to
    /// max(1, |coordinate|), in internal coordinates
    pub x_tol: f64,
    /// Spread of function values, relative to max(1, |best value|)
    pub f_tol: f64,
    /// Dimension-dependent coefficients (Gao & Han), recommended for many parameters
    pub adaptive: bool,
    /// Relative size of the first simplex step (0.05 = 5%)
    pub initial_step: f64,
}

impl Default for NelderMead {
    fn default() -> Self {
        Self {
            max_iter: 10_000,
            max_evals: 20_000,
            x_tol: 1e-4,
            f_tol: 1e-6,
            adaptive: true,
            initial_step: 0.05,
        }
    }
}

struct Run {
    u: Vec<f64>,
    f: f64,
    evals: usize,
    iters: usize,
    stop: StopReason,
}

/// c + t * (p - c)
fn lerp(c: &[f64], p: &[f64], t: f64) -> Vec<f64> {
    c.iter().zip(p).map(|(&c, &p)| c + t * (p - c)).collect()
}

impl NelderMead {
    /// Minimise `f` starting from `x0` (external values, must be inside `bounds`).
    ///
    /// `on_progress(iter, evals, best_x, best_f)` is called at the start of every
    /// iteration (and once more when the search stops); use it for checkpoints.
    /// Non-finite objective values are treated as `f64::MAX`.
    pub fn minimize<F, P>(
        &self,
        mut f: F,
        x0: &[f64],
        bounds: &Bounds,
        mut on_progress: P,
    ) -> Result<Solution>
    where
        F: FnMut(&[f64]) -> f64,
        P: FnMut(usize, usize, &[f64], f64),
    {
        ensure!(!x0.is_empty(), "no parameters to optimise");
        ensure!(x0.len() == bounds.len(), "start has {} values but there are {} bounds", x0.len(), bounds.len());
        ensure!(self.initial_step > 0.0, "initial_step must be positive");
        ensure!(self.max_iter > 0 && self.max_evals > 0, "iteration and evaluation limits must be positive");

        let u0 = bounds.to_internal(x0)?;
        let run = self.run(
            |u: &[f64]| {
                let v = f(&bounds.to_external(u));
                if v.is_finite() { v } else { f64::MAX }
            },
            u0,
            |iter, evals, u, fv| on_progress(iter, evals, &bounds.to_external(u), fv),
        );

        Ok(Solution {
            x: bounds.to_external(&run.u),
            f: run.f,
            evals: run.evals,
            iters: run.iters,
            stop: run.stop,
        })
    }

    fn converged(&self, simplex: &[(Vec<f64>, f64)]) -> bool {
        let best_u = &simplex[0].0;
        let best_f = simplex[0].1;

        let mut f_spread = 0.0f64;
        let mut x_spread = 0.0f64;
        for (v, fv) in simplex {
            f_spread = f_spread.max((fv - best_f).abs());
            for (a, b) in v.iter().zip(best_u.iter()) {
                x_spread = x_spread.max((a - b).abs() / b.abs().max(1.0));
            }
        }
        f_spread <= self.f_tol * best_f.abs().max(1.0) && x_spread <= self.x_tol
    }

    /// Nelder-Mead in internal (unbounded) coordinates
    fn run<E, P>(&self, mut eval: E, u0: Vec<f64>, mut progress: P) -> Run
    where
        E: FnMut(&[f64]) -> f64,
        P: FnMut(usize, usize, &[f64], f64),
    {
        let n = u0.len();
        let nf = n as f64;
        // reflection, expansion, contraction, shrink
        let (alpha, beta, gamma, delta) = if self.adaptive && n >= 2 {
            (1.0, 1.0 + 2.0 / nf, 0.75 - 1.0 / (2.0 * nf), 1.0 - 1.0 / nf)
        } else {
            (1.0, 2.0, 0.5, 0.5)
        };

        // initial simplex: x0 plus one vertex per coordinate
        let mut simplex: Vec<(Vec<f64>, f64)> = Vec::with_capacity(n + 1);
        let f0 = eval(&u0);
        simplex.push((u0.clone(), f0));
        for i in 0..n {
            let mut v = u0.clone();
            v[i] = if v[i] != 0.0 { v[i] * (1.0 + self.initial_step) } else { ZERO_STEP };
            let fv = eval(&v);
            simplex.push((v, fv));
        }
        let mut evals = n + 1;
        let mut iters = 0usize;

        let stop = loop {
            simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
            progress(iters, evals, &simplex[0].0, simplex[0].1);

            if self.converged(&simplex) {
                break StopReason::Converged;
            }
            if iters >= self.max_iter {
                break StopReason::MaxIter;
            }
            if evals >= self.max_evals {
                break StopReason::MaxEvals;
            }

            // centroid of all vertices except the worst
            let mut c = vec![0.0; n];
            for (v, _) in &simplex[..n] {
                for i in 0..n {
                    c[i] += v[i] / nf;
                }
            }
            let worst = simplex[n].clone();
            let f_best = simplex[0].1;
            let f_second = simplex[n - 1].1;

            let xr = lerp(&c, &worst.0, -alpha); // reflection
            let fr = eval(&xr);
            evals += 1;

            let mut shrink = false;
            if fr < f_best {
                let xe = lerp(&c, &xr, beta); // expansion
                let fe = eval(&xe);
                evals += 1;
                simplex[n] = if fe < fr { (xe, fe) } else { (xr, fr) };
            } else if fr < f_second {
                simplex[n] = (xr, fr);
            } else {
                let outside = fr < worst.1;
                let xc = if outside {
                    lerp(&c, &xr, gamma) // outside contraction
                } else {
                    lerp(&c, &worst.0, gamma) // inside contraction
                };
                let fc = eval(&xc);
                evals += 1;
                let accept = if outside { fc <= fr } else { fc < worst.1 };
                if accept {
                    simplex[n] = (xc, fc);
                } else {
                    shrink = true;
                }
            }

            if shrink {
                let best = simplex[0].0.clone();
                for j in 1..=n {
                    let v = lerp(&best, &simplex[j].0, delta);
                    let fv = eval(&v);
                    evals += 1;
                    simplex[j] = (v, fv);
                }
            }
            iters += 1;
        };

        Run { u: simplex[0].0.clone(), f: simplex[0].1, evals, iters, stop }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tight() -> NelderMead {
        NelderMead {
            max_iter: 100_000,
            max_evals: 100_000,
            x_tol: 1e-9,
            f_tol: 1e-14,
            adaptive: true,
            initial_step: 0.05,
        }
    }

    fn rosenbrock(x: &[f64]) -> f64 {
        let (a, b) = (x[0], x[1]);
        (1.0 - a).powi(2) + 100.0 * (b - a * a).powi(2)
    }

    #[test]
    fn transforms_round_trip() {
        let bounds = Bounds::new(vec![
            (0.0, 1.0),
            (0.0, f64::INFINITY),
            (f64::NEG_INFINITY, f64::INFINITY),
            (f64::NEG_INFINITY, 5.0),
            (2.0, 7.0),
        ])
        .unwrap();
        let x = [0.3, 4.0, -2.5, 1.5, 6.0];
        let back = bounds.to_external(&bounds.to_internal(&x).unwrap());
        for (a, b) in x.iter().zip(&back) {
            assert!((a - b).abs() < 1e-9, "{a} vs {b}");
        }

        // any internal vector lands inside the bounds
        let wild = bounds.to_external(&[10.0, -3.0, 7.0, 2.0, -100.0]);
        assert!(wild[0] >= 0.0 && wild[0] <= 1.0);
        assert!(wild[1] >= 0.0);
        assert!(wild[3] <= 5.0);
        assert!(wild[4] >= 2.0 && wild[4] <= 7.0);
    }

    #[test]
    fn bad_bounds_and_points_are_errors() {
        assert!(Bounds::new(vec![(1.0, 1.0)]).is_err());
        assert!(Bounds::new(vec![(2.0, 1.0)]).is_err());
        let b = Bounds::new(vec![(0.0, 1.0)]).unwrap();
        assert!(b.to_internal(&[2.0]).is_err());
        assert!(b.to_internal(&[0.5, 0.5]).is_err());
        let nm = NelderMead::default();
        assert!(nm.minimize(|_| 0.0, &[0.5, 0.5], &b, |_, _, _, _| {}).is_err());
    }

    #[test]
    fn quadratic_5d() {
        let target = [1.0, 2.0, 3.0, 4.0, 5.0];
        let f = |x: &[f64]| x.iter().zip(target.iter()).map(|(a, b)| (a - b) * (a - b)).sum::<f64>();
        let sol = tight().minimize(f, &[1.5; 5], &Bounds::unbounded(5), |_, _, _, _| {}).unwrap();
        println!("5-d quadratic: f = {:e}, evals = {}, stop = {:?}", sol.f, sol.evals, sol.stop);
        for (a, b) in sol.x.iter().zip(target.iter()) {
            assert!((a - b).abs() < 1e-3, "{a} vs {b}");
        }
    }

    #[test]
    fn rosenbrock_2d() {
        let nm = NelderMead { adaptive: false, ..tight() };
        let sol = nm.minimize(rosenbrock, &[-1.2, 1.0], &Bounds::unbounded(2), |_, _, _, _| {}).unwrap();
        println!("Rosenbrock: x = {:?}, f = {:e}, evals = {}", sol.x, sol.f, sol.evals);
        assert!((sol.x[0] - 1.0).abs() < 1e-3 && (sol.x[1] - 1.0).abs() < 1e-3);
    }

    #[test]
    fn default_settings_report_convergence() {
        let f = |x: &[f64]| (x[0] - 1.0).powi(2) + (x[1] - 2.0).powi(2);
        let sol = NelderMead::default()
            .minimize(f, &[0.5, 0.5], &Bounds::unbounded(2), |_, _, _, _| {})
            .unwrap();
        assert_eq!(sol.stop, StopReason::Converged);
    }

    #[test]
    fn max_evals_stops() {
        let nm = NelderMead { max_evals: 50, adaptive: false, ..tight() };
        let sol = nm.minimize(rosenbrock, &[-1.2, 1.0], &Bounds::unbounded(2), |_, _, _, _| {}).unwrap();
        assert_eq!(sol.stop, StopReason::MaxEvals);
        assert!(sol.evals >= 50 && sol.evals <= 60, "evals = {}", sol.evals);
    }

    #[test]
    fn two_sided_bound_is_respected() {
        // true minimum at x = 5, but x must stay in [0, 2]
        let mut max_seen = f64::NEG_INFINITY;
        let f = |x: &[f64]| {
            max_seen = max_seen.max(x[0]);
            (x[0] - 5.0).powi(2)
        };
        let bounds = Bounds::new(vec![(0.0, 2.0)]).unwrap();
        let sol = tight().minimize(f, &[1.0], &bounds, |_, _, _, _| {}).unwrap();
        println!("two-sided bound: x = {:?}", sol.x);
        assert!((sol.x[0] - 2.0).abs() < 1e-3);
        assert!(max_seen <= 2.0);
    }

    #[test]
    fn lower_bound_only() {
        // true minimum at x = -3, but x must stay >= 0
        let f = |x: &[f64]| (x[0] + 3.0).powi(2);
        let bounds = Bounds::new(vec![(0.0, f64::INFINITY)]).unwrap();
        let sol = tight().minimize(f, &[1.0], &bounds, |_, _, _, _| {}).unwrap();
        println!("lower bound: x = {:?}", sol.x);
        assert!(sol.x[0] >= 0.0 && sol.x[0] < 1e-3);
    }

    #[test]
    fn same_input_same_result() {
        let run = || {
            tight()
                .minimize(rosenbrock, &[-1.2, 1.0], &Bounds::unbounded(2), |_, _, _, _| {})
                .unwrap()
        };
        let (a, b) = (run(), run());
        assert_eq!(a.x, b.x);
        assert_eq!(a.f, b.f);
        assert_eq!(a.evals, b.evals);
    }
}
