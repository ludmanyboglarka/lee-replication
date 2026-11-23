//! Bounded Nelder-Mead (hand-written).
//!
//! * `Bounds` maps between the real parameter space (with bounds) and an internal,
//!   unbounded space, the same trick as MATLAB's `fminsearchbnd`:
//!     two-sided [lo, hi]: x = lo + (sin(u) + 1) / 2 * (hi - lo)
//!     lower bound only  : x = lo + u^2
//!     upper bound only  : x = hi - u^2
//!     no bounds         : x = u
//! * `nelder_mead` minimises a function over the unbounded space.
//! * `nelder_mead_bounded` does the mapping for you: the start point, the function
//!   arguments, the progress callback and the result are all in the real space.
//!
//! Do not start exactly on a bound: the start is then at u = 0, where the initial
//! simplex step falls back to a tiny absolute size.

use anyhow::{ensure, Result};

// ---------------------------------------------------------------------------
// Bounds
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct Bounds {
    lower: Vec<f64>,
    upper: Vec<f64>,
}

impl Bounds {
    /// One (lower, upper) pair per parameter; use +/- infinity for "no bound".
    pub fn new(pairs: &[(f64, f64)]) -> Result<Self> {
        let mut lower = Vec::with_capacity(pairs.len());
        let mut upper = Vec::with_capacity(pairs.len());
        for (i, &(lo, hi)) in pairs.iter().enumerate() {
            ensure!(
                !lo.is_nan() && !hi.is_nan() && lo < hi,
                "bounds of parameter {i} must satisfy lower < upper, got ({lo}, {hi})"
            );
            lower.push(lo);
            upper.push(hi);
        }
        Ok(Bounds { lower, upper })
    }

    pub fn len(&self) -> usize {
        self.lower.len()
    }

    pub fn is_empty(&self) -> bool {
        self.lower.is_empty()
    }

    /// Real values -> internal values. Errors if a value is outside its bounds.
    pub fn to_internal(&self, x: &[f64]) -> Result<Vec<f64>> {
        ensure!(
            x.len() == self.len(),
            "expected {} parameters, got {}",
            self.len(),
            x.len()
        );
        let mut u = Vec::with_capacity(x.len());
        for (i, &xi) in x.iter().enumerate() {
            let (lo, hi) = (self.lower[i], self.upper[i]);
            ensure!(
                xi >= lo && xi <= hi,
                "parameter {i} = {xi} is outside its bounds [{lo}, {hi}]"
            );
            u.push(match (lo.is_finite(), hi.is_finite()) {
                (true, true) => (2.0 * (xi - lo) / (hi - lo) - 1.0).clamp(-1.0, 1.0).asin(),
                (true, false) => (xi - lo).sqrt(),
                (false, true) => (hi - xi).sqrt(),
                (false, false) => xi,
            });
        }
        Ok(u)
    }

    /// Internal values -> real values. Always inside the bounds.
    pub fn to_external(&self, u: &[f64]) -> Vec<f64> {
        u.iter()
            .enumerate()
            .map(|(i, &ui)| {
                let (lo, hi) = (self.lower[i], self.upper[i]);
                match (lo.is_finite(), hi.is_finite()) {
                    (true, true) => (lo + (ui.sin() + 1.0) / 2.0 * (hi - lo)).clamp(lo, hi),
                    (true, false) => lo + ui * ui,
                    (false, true) => hi - ui * ui,
                    (false, false) => ui,
                }
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Nelder-Mead
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub struct NelderMeadSettings {
    pub max_iter: usize,
    pub max_evals: usize,
    /// Converged when every vertex is within `x_tol * max(1, |best|)` of the best, per coordinate...
    pub x_tol: f64,
    /// ... and every vertex value is within `f_tol` of the best value.
    pub f_tol: f64,
    /// Coefficients that depend on the number of parameters (better in high dimensions).
    pub adaptive: bool,
    /// Relative size of the initial simplex steps (0.05 = 5% of each value).
    pub initial_step: f64,
    /// Step used for a coordinate that is exactly zero.
    pub zero_step: f64,
}

impl Default for NelderMeadSettings {
    fn default() -> Self {
        NelderMeadSettings {
            max_iter: 100_000,
            max_evals: 20_000,
            x_tol: 1e-4,
            f_tol: 1e-4,
            adaptive: true,
            initial_step: 0.05,
            zero_step: 0.00025,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopReason {
    Converged,
    MaxIter,
    MaxEvals,
}

#[derive(Debug, Clone)]
pub struct NmResult {
    pub x: Vec<f64>,
    pub fval: f64,
    pub n_evals: usize,
    pub n_iters: usize,
    pub stop: StopReason,
}

/// Passed to the callback after every iteration (use it to print or checkpoint).
#[derive(Debug)]
pub struct Progress<'a> {
    pub iter: usize,
    pub evals: usize,
    /// Best value so far.
    pub fval: f64,
    /// Best point so far.
    pub x: &'a [f64],
}

/// Counts evaluations and turns NaN into +infinity so values can always be ordered.
struct Evaluator<F> {
    f: F,
    count: usize,
}

impl<F: FnMut(&[f64]) -> f64> Evaluator<F> {
    fn eval(&mut self, x: &[f64]) -> f64 {
        self.count += 1;
        let v = (self.f)(x);
        if v.is_nan() {
            f64::INFINITY
        } else {
            v
        }
    }
}

/// c + t * (p - c)
fn point(c: &[f64], p: &[f64], t: f64) -> Vec<f64> {
    c.iter().zip(p.iter()).map(|(&ci, &pi)| ci + t * (pi - ci)).collect()
}

fn converged(verts: &[(Vec<f64>, f64)], x_tol: f64, f_tol: f64) -> bool {
    let best = &verts[0];
    verts.iter().skip(1).all(|v| {
        (v.1 - best.1).abs() <= f_tol
            && v.0
                .iter()
                .zip(best.0.iter())
                .all(|(a, b)| (a - b).abs() <= x_tol * b.abs().max(1.0))
    })
}

/// Minimise `f` over the whole real space, starting at `x0`.
pub fn nelder_mead<F, C>(
    f: F,
    x0: &[f64],
    s: &NelderMeadSettings,
    mut callback: C,
) -> Result<NmResult>
where
    F: FnMut(&[f64]) -> f64,
    C: FnMut(&Progress<'_>),
{
    let n = x0.len();
    ensure!(n > 0, "need at least one parameter");
    ensure!(x0.iter().all(|v| v.is_finite()), "start point must be finite");

    let (alpha, gamma, rho, sigma) = if s.adaptive && n >= 2 {
        let nf = n as f64;
        (1.0, 1.0 + 2.0 / nf, 0.75 - 1.0 / (2.0 * nf), 1.0 - 1.0 / nf)
    } else {
        (1.0, 2.0, 0.5, 0.5)
    };

    let mut ev = Evaluator { f, count: 0 };

    // initial simplex: x0 plus one nudged copy per coordinate
    let mut verts: Vec<(Vec<f64>, f64)> = Vec::with_capacity(n + 1);
    let f0 = ev.eval(x0);
    verts.push((x0.to_vec(), f0));
    for i in 0..n {
        let mut v = x0.to_vec();
        v[i] = if v[i] != 0.0 { v[i] * (1.0 + s.initial_step) } else { s.zero_step };
        let fv = ev.eval(&v);
        verts.push((v, fv));
    }

    let mut iters = 0usize;
    let stop;
    loop {
        verts.sort_by(|a, b| a.1.total_cmp(&b.1));
        callback(&Progress { iter: iters, evals: ev.count, fval: verts[0].1, x: &verts[0].0 });

        if converged(&verts, s.x_tol, s.f_tol) {
            stop = StopReason::Converged;
            break;
        }
        if iters >= s.max_iter {
            stop = StopReason::MaxIter;
            break;
        }
        if ev.count >= s.max_evals {
            stop = StopReason::MaxEvals;
            break;
        }
        iters += 1;

        // centroid of all vertices except the worst
        let mut c = vec![0.0; n];
        for v in &verts[..n] {
            for j in 0..n {
                c[j] += v.0[j];
            }
        }
        for cj in c.iter_mut() {
            *cj /= n as f64;
        }

        let f_best = verts[0].1;
        let f_second = verts[n - 1].1;
        let f_worst = verts[n].1;

        // reflect the worst vertex through the centroid
        let xr = point(&c, &verts[n].0, -alpha);
        let fr = ev.eval(&xr);

        if fr < f_best {
            // expand
            let xe = point(&c, &xr, gamma);
            let fe = ev.eval(&xe);
            verts[n] = if fe < fr { (xe, fe) } else { (xr, fr) };
        } else if fr < f_second {
            verts[n] = (xr, fr);
        } else {
            let shrink;
            if fr < f_worst {
                // outside contraction
                let xc = point(&c, &xr, rho);
                let fc = ev.eval(&xc);
                if fc <= fr {
                    verts[n] = (xc, fc);
                    shrink = false;
                } else {
                    shrink = true;
                }
            } else {
                // inside contraction
                let xc = point(&c, &verts[n].0, rho);
                let fc = ev.eval(&xc);
                if fc < f_worst {
                    verts[n] = (xc, fc);
                    shrink = false;
                } else {
                    shrink = true;
                }
            }
            if shrink {
                let best = verts[0].0.clone();
                for v in verts.iter_mut().skip(1) {
                    for j in 0..n {
                        v.0[j] = best[j] + sigma * (v.0[j] - best[j]);
                    }
                    let fv = ev.eval(&v.0);
                    v.1 = fv;
                }
            }
        }
    }

    verts.sort_by(|a, b| a.1.total_cmp(&b.1));
    Ok(NmResult {
        x: verts[0].0.clone(),
        fval: verts[0].1,
        n_evals: ev.count,
        n_iters: iters,
        stop,
    })
}

/// Minimise `f` subject to `bounds`. Start point, function arguments, callback and
/// result are all in the real (bounded) parameter space.
pub fn nelder_mead_bounded<F, C>(
    mut f: F,
    x0: &[f64],
    bounds: &Bounds,
    s: &NelderMeadSettings,
    mut callback: C,
) -> Result<NmResult>
where
    F: FnMut(&[f64]) -> f64,
    C: FnMut(&Progress<'_>),
{
    ensure!(
        x0.len() == bounds.len(),
        "start point has {} values but there are {} bounds",
        x0.len(),
        bounds.len()
    );
    let u0 = bounds.to_internal(x0)?;

    let mut res = nelder_mead(
        |u: &[f64]| f(&bounds.to_external(u)),
        &u0,
        s,
        |p: &Progress<'_>| {
            let ext = bounds.to_external(p.x);
            callback(&Progress { iter: p.iter, evals: p.evals, fval: p.fval, x: &ext });
        },
    )?;
    res.x = bounds.to_external(&res.x);
    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tight() -> NelderMeadSettings {
        NelderMeadSettings {
            max_iter: 200_000,
            max_evals: 200_000,
            x_tol: 1e-9,
            f_tol: 1e-14,
            ..Default::default()
        }
    }

    fn mixed_bounds() -> Bounds {
        Bounds::new(&[
            (0.0, 1.0),
            (0.0, f64::INFINITY),
            (f64::NEG_INFINITY, 5.0),
            (f64::NEG_INFINITY, f64::INFINITY),
        ])
        .unwrap()
    }

    #[test]
    fn bounds_round_trip() {
        let b = mixed_bounds();
        let x = [0.3, 2.0, 1.0, -7.0];
        let u = b.to_internal(&x).unwrap();
        let back = b.to_external(&u);
        for (a, c) in x.iter().zip(back.iter()) {
            assert!((a - c).abs() < 1e-9, "{a} vs {c}");
        }
    }

    #[test]
    fn any_internal_value_maps_inside_the_bounds() {
        let b = mixed_bounds();
        for &u in &[-100.0, -10.0, -3.0, -0.5, 0.0, 1.0, 2.5, 10.0, 100.0] {
            let x = b.to_external(&[u, u, u, u]);
            assert!(x[0] >= 0.0 && x[0] <= 1.0);
            assert!(x[1] >= 0.0);
            assert!(x[2] <= 5.0);
        }
    }

    #[test]
    fn bad_bounds_and_outside_values_are_errors() {
        assert!(Bounds::new(&[(1.0, 1.0)]).is_err());
        assert!(Bounds::new(&[(2.0, 1.0)]).is_err());
        let b = mixed_bounds();
        assert!(b.to_internal(&[1.5, 2.0, 1.0, 0.0]).is_err());
        assert!(b.to_internal(&[0.5, -1.0, 1.0, 0.0]).is_err());
        assert!(b.to_internal(&[0.5, 1.0, 6.0, 0.0]).is_err());
        assert!(b.to_internal(&[0.5, 1.0]).is_err());
    }

    #[test]
    fn finds_minimum_of_a_quadratic() {
        let target = [1.0, -2.0, 3.0, 0.5, 4.0];
        let f = |x: &[f64]| x.iter().zip(target.iter()).map(|(a, t)| (a - t).powi(2)).sum::<f64>();
        for adaptive in [true, false] {
            let s = NelderMeadSettings { adaptive, ..tight() };
            let r = nelder_mead(f, &[1.0; 5], &s, |_| {}).unwrap();
            assert!(r.fval < 1e-6, "adaptive={adaptive}: fval {}", r.fval);
            for (a, t) in r.x.iter().zip(target.iter()) {
                assert!((a - t).abs() < 1e-3);
            }
        }
    }

    #[test]
    fn finds_minimum_of_rosenbrock() {
        let f = |x: &[f64]| 100.0 * (x[1] - x[0] * x[0]).powi(2) + (1.0 - x[0]).powi(2);
        let r = nelder_mead(f, &[-1.2, 1.0], &tight(), |_| {}).unwrap();
        assert!((r.x[0] - 1.0).abs() < 1e-3 && (r.x[1] - 1.0).abs() < 1e-3, "{:?}", r.x);
        assert_eq!(r.stop, StopReason::Converged);
    }

    #[test]
    fn respects_bounds_when_the_optimum_is_outside() {
        // unconstrained minimum at (3, -2); allowed region x in [0, 1], y in [-1, 5]
        let bounds = Bounds::new(&[(0.0, 1.0), (-1.0, 5.0)]).unwrap();
        let mut all_inside = true;
        let r = nelder_mead_bounded(
            |x: &[f64]| {
                if x[0] < 0.0 || x[0] > 1.0 || x[1] < -1.0 || x[1] > 5.0 {
                    all_inside = false;
                }
                (x[0] - 3.0).powi(2) + (x[1] + 2.0).powi(2)
            },
            &[0.5, 0.0],
            &bounds,
            &tight(),
            |_| {},
        )
        .unwrap();
        assert!(all_inside);
        assert!((r.x[0] - 1.0).abs() < 1e-3, "{:?}", r.x);
        assert!((r.x[1] + 1.0).abs() < 1e-3, "{:?}", r.x);
    }

    #[test]
    fn stops_at_the_evaluation_limit() {
        let f = |x: &[f64]| 100.0 * (x[1] - x[0] * x[0]).powi(2) + (1.0 - x[0]).powi(2);
        let s = NelderMeadSettings { max_evals: 10, ..tight() };
        let r = nelder_mead(f, &[-1.2, 1.0], &s, |_| {}).unwrap();
        assert_eq!(r.stop, StopReason::MaxEvals);
    }

    #[test]
    fn callback_sees_a_never_worsening_best_value() {
        let f = |x: &[f64]| (x[0] - 2.0).powi(2) + (x[1] + 1.0).powi(2);
        let mut last = f64::INFINITY;
        let mut calls = 0usize;
        nelder_mead(f, &[3.0, 3.0], &tight(), |p: &Progress<'_>| {
            assert!(p.fval <= last);
            last = p.fval;
            calls += 1;
        })
        .unwrap();
        assert!(calls > 5);
    }

    #[test]
    fn rejects_empty_or_non_finite_start() {
        let f = |_: &[f64]| 0.0;
        assert!(nelder_mead(f, &[], &tight(), |_| {}).is_err());
        assert!(nelder_mead(f, &[f64::NAN], &tight(), |_| {}).is_err());
    }
}
