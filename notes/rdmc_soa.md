## What the function does
### Inputs
`A`: how much weight the automatic process has at the start of the trial.
`k`: how fast that weight decays.
`v_a`: the automatic drift. The caller passes it positive for congruent trials and negative for incongruent ones.
`v_c`: the controlled drift, pointing to the correct response.
`c`: the starting point (bias), 0 in your models.
- Fixed constants. The noise scale σ is 0.1. Time runs in 1 ms steps for 1.5 s, so 1,500 steps, with time measured in seconds.

### Each step.
    1. Weights. At time t, the automatic weight is A·e^(−k·t) and the controlled weight is 1 minus that.
    2. Drift. The total drift is the weighted mix of the two: w_a·v_a + w_c·v_c. Early on it's dominated by the automatic drift, and later by the controlled drift.
    3. Increment. Each step moves the trace by drift × dt plus Gaussian noise of size σ·√dt. This is the standard step for simulating a diffusion process.
    4. Position. The trace is the running sum of the increments plus the start point c.

### Output. 
- An array of 1,500 positions, where position i is the state at i ms. The function knows nothing about boundaries. It never says whether the trial ended or which response won.

### Dead code. 
- The lines that compute dX_a, X_a, dX_c and X_c are never used. 
- -> they are not translated to rust

### What changes in Rust
- The boundary logic lives elsewhere in MATLAB (in the objective file, after all traces are made). In Rust, you'll fold it in:

    - Step through time one step at a time, update the position, and check the boundaries.
    - Stop at the first crossing and return which boundary and at which step.
    - You never need to store the path, and trials end early, which is where the speedup comes from.
    - The non-decision time (Ter) is added later, outside this function.

- Keep the details exactly as the MATLAB: the first step uses t = 0.001, not 0, and drift and noise are per second.


# What `model.rs` does 
- `#[derive(Debug, Clone, Copy, PartialEq, Eq)]`: 
    - Debug: print the struct for inspection. e.g. `println!("{:?}", p)` shows TrialParams `{ a_weight: 0.5, k: 8.0, ... }`. Without it, printing a TrialParams is a compile error.
    - Clone: lets me make an explicit duplicate with p.clone()
    - Copy: makes the duplication automatic
    - PartialEq means you can ask whether two of them are identical.
    - Eq is a stricter version of PartialEq. PartialEq gives you ==. Eq adds a promise that the comparison is always well-behaved: every value is equal to itself, so x == x is always true.

- next step: 
pub enum Outcome {
    Correct { step: u32 }, 
    Error { step: u32 },
    Timeout,
}

    - An enum is a value that must be exactly one of several named options, and each option can carry data.
    - Correct { step: 312 } means the trial hit the upper boundary at 312 ms.
    - Error { step: 405 } means it hit the lower boundary at 405 ms.
    - Timeout means it hit neither boundary within 1500 ms. It carries no data.

### simulation function: simulates one trial: 
1. the signature: 
- `pub fn simulate_trial` declares a function that other files can use.
- `p: &TrialParams` is the parameter bundle. The & means "borrow it, don't take it". The function can look at p but can't change it, and the caller keeps ownership. This is why Copy wasn't needed here.
- `rng: &mut R` is the random number generator. &mut means a borrow that the function is allowed to change. Generating a random number updates the generator's internal state, so it has to be changeable. In R terms, the generator is the hidden `.Random.seed`, passed in explicitly instead of living globally.
- -> Outcome says what comes back: one of Correct, Error or Timeout.

# `data.rs`
- What derive(Deserialize) generates
The compiler reads your struct definition (field names and types) and writes the conversion code for you. In effect it produces the logic:
find the column named seq and read it as a String
find the column named soa and parse it as an i32
find e1 and parse it as an f64
and so on for every field
- Where serde and csv fit
The work is split between two crates:
serde defines the general idea of "turn data into a Rust struct", and provides the Deserialize derive.
csv knows the CSV format and uses that generated code in rdr.deserialize::<Condition>().

- Where it sits in the pipeline
R script  →  flanker_cse_observed.csv  →  data.rs  →  Vec<Condition>  →  objective.rs  (read + check)                (compare with simulation)
- data.rs is the boundary between the two languages. The CSV is the only contact point between R and Rust, and this file is the Rust side of that contact. If the CSV format changes, this is the only file that has to change.






