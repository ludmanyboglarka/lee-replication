Variant has two options, Simple (26 parameters) and Da (30). Each variant has n_params(), names() (like A_100, kc_0, a, std), bounds(), and expand().
One table per variant lists the parameter blocks in vector order. The names, bounds and positions all come from that table, so they can't drift apart.
expand returns one ConditionModel per condition, in the order of the data. For each condition, the SOA picks the position (100, 0, −50, −100 are 1 to 4) and the sequence label gives previous and current congruency. The current congruency picks kc or ki and the sign of the automatic drift. The previous congruency picks the da set, which only matters in the Da variant.
Bounds: A in [0, 1], the decay rates k ≥ 0, the boundary a ≥ 0.001, and Ter and std ≥ 0. The drift rates da and dc are unbounded.
The six tests
Sizes agree: the counts, names and bounds all line up, and the boundary, Ter and std land at the expected indices.
Wrong length: a vector of the wrong size returns an error.
Unknown labels: an unknown SOA or sequence label returns an error.
Hand-checked cases: one test for each variant. The parameter vector is 0, 1, 2, …, so every slot has a unique value. The test expands ci at SOA −50 for the simple model and ic at SOA 0 for the Da model, then checks that every field came from the right slot.
Nesting: a Da vector whose two da sets are equal gives exactly the same 16 models as the matching simple vector.
The file uses panic! in one place, if a block name is missing from the table. That can only happen from a typo in the code itself, which the tests would catch, so I didn't treat it as a runtime error.
Next is optim.rs: a bounded Nelder-Mead that calls expand and deviance, with the coarse-then-fine trial counts and a way to save its progress. Want me to explain how that will work first?

