//! Model variants: parameter vector -> one `ConditionModel` per condition, plus names and bounds.
//!
//! Vector layouts (from the MATLAB wrapper), SOA order 100, 0, -50, -100:
//!
//!   Simple (26): A[4] kc[4] ki[4] da[4] dc[4] a Ter[4] std
//!   Da     (30): A[4] kc[4] ki[4] da_pc[4] da_pi[4] dc[4] a Ter[4] std
//!
//! Sequence labels: first letter = PREVIOUS trial, second = CURRENT trial.
//!   kc / ki : decay rate when the current trial is congruent / incongruent
//!   da      : automatic drift, + on congruent and - on incongruent current trials
//!   da_pc / da_pi (Da only): da after a congruent / incongruent previous trial

use anyhow::{anyhow, bail, ensure, Result};

use crate::data::Condition;
use crate::model::TrialParams;
use crate::objective::ConditionModel;

/// SOA of parameter position 1..4.
pub const SOAS: [i32; 4] = [100, 0, -50, -100];

/// Starting point of the evidence (the CSE objective uses c = 0).
const START: f64 = 0.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Variant {
    /// CSE model 1: one `da` set, 26 parameters.
    Simple,
    /// CSE model 2: `da` depends on the previous trial, 30 parameters.
    Da,
}

/// Previous and current congruency of a sequence label: (prev_congruent, cur_congruent).
fn parse_seq(seq: &str) -> Result<(bool, bool)> {
    match seq {
        "cc" => Ok((true, true)),
        "ci" => Ok((true, false)),
        "ic" => Ok((false, true)),
        "ii" => Ok((false, false)),
        _ => bail!("unknown sequence label {seq:?} (expected cc, ci, ic or ii)"),
    }
}

impl Variant {
    pub fn n_params(self) -> usize {
        match self {
            Variant::Simple => 26,
            Variant::Da => 30,
        }
    }

    /// Parameter blocks in vector order: (name, number of values).
    /// A block with 4 values has one value per SOA.
    fn blocks(self) -> Vec<(&'static str, usize)> {
        match self {
            Variant::Simple => vec![
                ("A", 4), ("kc", 4), ("ki", 4), ("da", 4), ("dc", 4),
                ("a", 1), ("Ter", 4), ("std", 1),
            ],
            Variant::Da => vec![
                ("A", 4), ("kc", 4), ("ki", 4), ("da_pc", 4), ("da_pi", 4), ("dc", 4),
                ("a", 1), ("Ter", 4), ("std", 1),
            ],
        }
    }

    /// Index in the vector where a block starts.
    fn start_of(self, block: &str) -> usize {
        let mut pos = 0;
        for (name, count) in self.blocks() {
            if name == block {
                return pos;
            }
            pos += count;
        }
        panic!("block {block:?} missing from the {self:?} layout");
    }

    /// One name per parameter, e.g. "A_100", "kc_0", "a", "std".
    pub fn names(self) -> Vec<String> {
        let mut out = Vec::new();
        for (name, count) in self.blocks() {
            if count == SOAS.len() {
                for soa in SOAS {
                    out.push(format!("{name}_{soa}"));
                }
            } else {
                out.push(name.to_string());
            }
        }
        out
    }

    /// (lower, upper) bound for every parameter; infinite where unbounded.
    pub fn bounds(self) -> Vec<(f64, f64)> {
        let mut out = Vec::new();
        for (name, count) in self.blocks() {
            let b = match name {
                "A" => (0.0, 1.0),
                "kc" | "ki" => (0.0, f64::INFINITY),
                "a" => (1e-3, f64::INFINITY),
                "Ter" | "std" => (0.0, f64::INFINITY),
                _ => (f64::NEG_INFINITY, f64::INFINITY), // da, da_pc, da_pi, dc
            };
            out.extend(std::iter::repeat(b).take(count));
        }
        out
    }

    /// Build the model of every condition, in the same order as `conditions`.
    pub fn expand(self, w: &[f64], conditions: &[Condition]) -> Result<Vec<ConditionModel>> {
        ensure!(
            w.len() == self.n_params(),
            "{:?} model needs {} parameters, got {}",
            self,
            self.n_params(),
            w.len()
        );

        let a_weight = self.start_of("A");
        let kc = self.start_of("kc");
        let ki = self.start_of("ki");
        let (da_pc, da_pi) = match self {
            Variant::Simple => {
                let d = self.start_of("da");
                (d, d)
            }
            Variant::Da => (self.start_of("da_pc"), self.start_of("da_pi")),
        };
        let dc = self.start_of("dc");
        let boundary = self.start_of("a");
        let ter = self.start_of("Ter");
        let std_idx = self.start_of("std");

        conditions
            .iter()
            .map(|c| -> Result<ConditionModel> {
                let m = SOAS
                    .iter()
                    .position(|&s| s == c.soa)
                    .ok_or_else(|| anyhow!("unknown SOA {} in condition {}", c.soa, c.seq))?;
                let (prev_cong, cur_cong) = parse_seq(&c.seq)?;

                let k = if cur_cong { w[kc + m] } else { w[ki + m] };
                let da = if prev_cong { w[da_pc + m] } else { w[da_pi + m] };
                let v_a = if cur_cong { da } else { -da };

                Ok(ConditionModel {
                    trial: TrialParams {
                        a_weight: w[a_weight + m],
                        k,
                        v_a,
                        v_c: w[dc + m],
                        start: START,
                        boundary: w[boundary],
                    },
                    ter: w[ter + m],
                    std: w[std_idx],
                })
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 16 conditions in CSV order: cc, ci, ic, ii, each with SOA 100, 0, -50, -100.
    fn all_conditions() -> Vec<Condition> {
        let mut v = Vec::new();
        for seq in ["cc", "ci", "ic", "ii"] {
            for soa in SOAS {
                v.push(Condition {
                    seq: seq.to_string(),
                    soa,
                    e1: 1.0, e2: 2.0, e3: 3.0, e4: 4.0, e5: 5.0, e6: 6.0,
                    f1: 1.0, f2: 1.0, f3: 1.0, f4: 1.0, f5: 1.0, f6: 1.0,
                    f_err: 1.0,
                });
            }
        }
        v
    }

    /// Parameter vector whose i-th value is i, so every slot is distinguishable.
    fn counting(n: usize) -> Vec<f64> {
        (0..n).map(|i| i as f64).collect()
    }

    #[test]
    fn sizes_names_and_bounds_agree() {
        for v in [Variant::Simple, Variant::Da] {
            assert_eq!(v.names().len(), v.n_params());
            assert_eq!(v.bounds().len(), v.n_params());
        }
        assert_eq!(Variant::Simple.n_params(), 26);
        assert_eq!(Variant::Da.n_params(), 30);

        let names = Variant::Simple.names();
        assert_eq!(names[0], "A_100");
        assert_eq!(names[20], "a");
        assert_eq!(names[21], "Ter_100");
        assert_eq!(names[25], "std");
        let names = Variant::Da.names();
        assert_eq!(names[16], "da_pi_100");
        assert_eq!(names[24], "a");
        assert_eq!(names[29], "std");
    }

    #[test]
    fn wrong_length_is_an_error() {
        let conds = all_conditions();
        assert!(Variant::Simple.expand(&counting(25), &conds).is_err());
        assert!(Variant::Da.expand(&counting(26), &conds).is_err());
    }

    #[test]
    fn unknown_labels_are_errors() {
        let mut conds = all_conditions();
        conds[0].soa = 42;
        assert!(Variant::Simple.expand(&counting(26), &conds).is_err());
        let mut conds = all_conditions();
        conds[0].seq = "xx".to_string();
        assert!(Variant::Simple.expand(&counting(26), &conds).is_err());
    }

    #[test]
    fn simple_hand_checked_case() {
        // "ci" at SOA -50: condition index 1*4 + 2 = 6, SOA position m = 2
        let models = Variant::Simple.expand(&counting(26), &all_conditions()).unwrap();
        assert_eq!(models.len(), 16);
        let m = models[6];
        assert_eq!(m.trial.a_weight, 2.0); // A[2]
        assert_eq!(m.trial.k, 10.0); // ki[2] (current incongruent)
        assert_eq!(m.trial.v_a, -14.0); // -da[2]
        assert_eq!(m.trial.v_c, 18.0); // dc[2]
        assert_eq!(m.trial.start, 0.0);
        assert_eq!(m.trial.boundary, 20.0); // a
        assert_eq!(m.ter, 23.0); // Ter[2]
        assert_eq!(m.std, 25.0);
    }

    #[test]
    fn da_hand_checked_case() {
        // "ic" at SOA 0: condition index 2*4 + 1 = 9, SOA position m = 1
        let models = Variant::Da.expand(&counting(30), &all_conditions()).unwrap();
        let m = models[9];
        assert_eq!(m.trial.a_weight, 1.0); // A[1]
        assert_eq!(m.trial.k, 5.0); // kc[1] (current congruent)
        assert_eq!(m.trial.v_a, 17.0); // +da_pi[1] (previous incongruent)
        assert_eq!(m.trial.v_c, 21.0); // dc[1]
        assert_eq!(m.trial.boundary, 24.0); // a
        assert_eq!(m.ter, 26.0); // Ter[1]
        assert_eq!(m.std, 29.0);
    }

    #[test]
    fn da_with_equal_da_sets_equals_simple() {
        let ws = counting(26);
        // Da vector: A kc ki | da, da | dc a Ter std
        let mut wd = Vec::new();
        wd.extend_from_slice(&ws[0..12]);
        wd.extend_from_slice(&ws[12..16]);
        wd.extend_from_slice(&ws[12..16]);
        wd.extend_from_slice(&ws[16..26]);
        assert_eq!(wd.len(), 30);

        let conds = all_conditions();
        let a = Variant::Simple.expand(&ws, &conds).unwrap();
        let b = Variant::Da.expand(&wd, &conds).unwrap();
        assert_eq!(format!("{a:?}"), format!("{b:?}"));
    }
}

