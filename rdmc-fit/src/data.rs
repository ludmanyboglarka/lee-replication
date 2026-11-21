//! Loader for `data/flanker_cse_observed.csv` (written by r/01_preprocess.R).
//!
//! One row per condition (sequence x SOA):
//!   seq, soa, e1..e6 (RT bin edges, ms), f1..f6 (observed correct-RT bin counts), f_err

use std::io::Read;
use std::path::Path;

use anyhow::{bail, Context, Result};
use serde::Deserialize;

#[derive(Debug, Clone, Deserialize)]
pub struct Condition {
    pub seq: String, // "cc", "ci", "ic", "ii" (first letter = previous trial)
    pub soa: i32,
    pub e1: f64,
    pub e2: f64,
    pub e3: f64,
    pub e4: f64,
    pub e5: f64,
    pub e6: f64,
    pub f1: f64,
    pub f2: f64,
    pub f3: f64,
    pub f4: f64,
    pub f5: f64,
    pub f6: f64,
    pub f_err: f64,
}

impl Condition {
    pub fn edges(&self) -> [f64; 6] {
        [self.e1, self.e2, self.e3, self.e4, self.e5, self.e6]
    }

    pub fn freqs(&self) -> [f64; 6] {
        [self.f1, self.f2, self.f3, self.f4, self.f5, self.f6]
    }

    /// Total observed frequency: six correct bins + error bin.
    pub fn total(&self) -> f64 {
        self.freqs().iter().sum::<f64>() + self.f_err
    }

    fn validate(&self) -> Result<()> {
        if !self.edges().windows(2).all(|w| w[0] < w[1]) {
            bail!("seq {} soa {}: bin edges are not increasing", self.seq, self.soa);
        }
        let mut all = self.freqs().to_vec();
        all.push(self.f_err);
        if all.iter().any(|f| !f.is_finite() || *f < 0.0) {
            bail!("seq {} soa {}: frequencies must be finite and >= 0", self.seq, self.soa);
        }
        Ok(())
    }
}

pub fn load(path: &Path) -> Result<Vec<Condition>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("cannot open observed-data file {}", path.display()))?;
    from_reader(file)
}

pub fn from_reader<R: Read>(reader: R) -> Result<Vec<Condition>> {
    let mut rdr = csv::Reader::from_reader(reader);
    let mut out = Vec::new();
    for (i, rec) in rdr.deserialize::<Condition>().enumerate() {
        let c = rec.with_context(|| format!("cannot parse data row {}", i + 1))?;
        c.validate().with_context(|| format!("invalid data row {}", i + 1))?;
        out.push(c);
    }
    if out.is_empty() {
        bail!("observed-data file has no rows");
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "seq,soa,e1,e2,e3,e4,e5,e6,f1,f2,f3,f4,f5,f6,f_err
cc,100,406.7,446.8,473.4,500.3,546.9,721.4,10.45,9.875,11.875,7.775,7.925,9.75,0.225
cc,0,407.8,442.4,465.8,493.0,543.1,723.4,9.7,12.525,6.55,12.475,10.375,8.125,0.3
";

    #[test]
    fn parses_sample() {
        let rows = from_reader(SAMPLE.as_bytes()).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].seq, "cc");
        assert_eq!(rows[1].soa, 0);
        assert_eq!(rows[0].edges()[5], 721.4);
        assert!((rows[0].total() - 57.875).abs() < 1e-9);
    }

    #[test]
    fn rejects_decreasing_edges() {
        let bad = SAMPLE.replace("406.7,446.8", "446.8,406.7");
        assert!(from_reader(bad.as_bytes()).is_err());
    }

    #[test]
    fn rejects_empty_file() {
        let header_only = SAMPLE.lines().next().unwrap();
        assert!(from_reader(header_only.as_bytes()).is_err());
    }
}
