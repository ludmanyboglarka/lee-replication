# What would the handwritten Nelder-Mead do? 
### The problem: 
- deviance(x) takes a vector of 26 (or 30) numbers and returns one number -> we want the vector that makes it the smallest 
- the method compares function values

### The simplex: 
- in n dimensions, NM keeps n+1 candidate points, this is a a simplex
-> in 2 dimensions: this is a triangle; and for us it is 27 points in 26 dimensions
- the algorithm moves the triangle across the landscape, away from bad regions and towards good ones
- Each point costs one full objective evaluation, which is 16 conditions of simulated trials

### The setup: 
- from starting vector -> NM builds the other n points by nudging one coordinate at a time
    (MATLAB's fminsearch nudges each coordinate by 5% of its own value, or by a tiny absolute amount (0.00025) if the value is zero)
    nudging a % keeps the step sensible even though our parameters live on very different scales (Ter is around 200, a is around 0.1)
- all n+1 points are evaluated and sorted from best to worst 

### One iteration
- Let the worst point be W. NM computes the centroid C of all the other points and tries to replace W. 
1. Reflect: 
    - try the mirror image of W through C, at C+(-W) -> "walk away from the worst point"
2. If the reflected point is the new best, expand. The direction is evidently good, so try going twice as far. Keep whichever of the two (reflected or expanded) is better
3. if its in the middle: better than the second-worst point, but not the best), just accept the reflection
4. if its worse than the second-worst but better than W, contract outside
    - Try a point halfway between C and the reflected point. Accept it if it improves on the reflection
5. If it is worse than W, contract inside
    - Pull every point except the best halfway towards the best one and re-evaluate them all. This is the expensive move, n evaluations in one go, and it means "this whole region is bad, zoom in"

textbook coefficients are 1 (reflect), 2 (expand), 0.5 (contract) and 0.5 (shrink)
    -> normal iteration: costs one or two evaluations, and a shrink costs n

### Weakness in high dimensions
- textbook versin works well in 2-5 dimensions -> degrades notably beyond about 10
- simplex flattens into a sliver and crawl (lol such a beautiful way of saying)

### Stopping: 
- NM has no gradient to tell when it has arrived, so it stops on one of these: 
    - the simplex has become tiny, meaning the points are all very close together
    - the function values across the simplex are almost identical
    - an evaluation or iteration budget is exhausted

- my evaluation: 
    - 26 parameters, 100 iterations and 3500 evaluations -> 130 evaluations per parameter 

- function-value criterion -> tolerance on parameters should be relative or scale-aware

### On our particular function: 
- fixed seed: `deviance(x)` is deterministic -> NM is not being misled by noise 
- The function does have tiny jumps, because a simulated trial can hop from one RT bin to another as a parameter changes
- With 8,000 or 50,000 trials those jumps are small. NM only compares values, so it copes with them much better than gradient methods would

### What it reports: 
- The best vector found, its deviance, the number of evaluations and iterations, and the reason it stopped
- also reports the best point ever seen 
- can call a progress function every so often -> checkpointing will be put in 
- vertices are evaluated one at a time, the parallelism lives inside the objective -- across the conditions

### Where it can fail
- premature collapse: the simplex shrinks to a degenerate shape and stops while still far from a true minimum
- Local minima: NM only goes downhill from where it starts 
- Flat directions: if several parameters trade off (for example A and k, or da and dc), the objective barely changes along that direction. NM wanders along the valley floor and stops with tiny improvements. This is often a property of the model, not a bug

# Handling the bounds: 
### The mappings: 
##### Two-sided [lo, hi]: 
- for example A in `[0, 1]: x = lo + (sin(u) + 1)/2 · (hi − lo)` 
- Any real u gives a valid x, so the optimizer can never leave the range. The inverse, to convert the start vector, uses the arcsine.y
##### Lower bound only 
- for example k ≥ 0, Ter ≥ 0: `x = lo + u²`
- Any u gives x ≥ lo. The original used this, so it keeps the behaviour faithful. The alternative lo + exp(u) is one-to-one, but it can never reach the bound and fails on a start exactly at the bound.
##### No bounds (da, dc): identity.

## Consequences: 
- the optimizer lives inside the internal space (5% step is 5% of the internal value)
- The sin and square maps are many-to-one (u and −u, or u and u + 2π, give the same x). That is harmless for NM, which only compares function values
- Unbounded-above parameters such as Ter, k and std can drift to absurd values. You may want finite upper limits (for example Ter below a few hundred ms), which can be set in the config as optional overrides.

### Tests we'd want:
- Converting to internal and back should return the original vector, any internal vector should map inside the bounds, and NM should find the minimum of a simple known function. A second test should use a function whose true minimum lies outside the bounds, and the result should end up at the bound

# Starting value and a config file 
- NM is local, and our landscape is hard 
    - It finds a low point near where it begins. With 26–30 parameters that trade off against each other, the deviance surface has valleys and possibly several local minima, so where you start affects where you end up
- Some starts are simply useless
- The o.g. code had its own starts 

### Why a config file and not hard-coded numbers: 
    - Reproducibility. The starting point is part of the definition of the analysis. If it lives in configs/ and is committed, anyone can see exactly what produced a result.
    - No recompiling. You can try a different start by editing a text file.
    - Different models need different starts. The "da" model is nested in the simple one. If you start it from the simple model's best solution (with both da sets equal), it begins at exactly that deviance and can only improve. This also makes the deviance difference between the two models meaningful for the chi-square test. A config file makes that handoff easy.
    - It collects the other run settings in one place. Trial counts for each stage, the seed, the number of restarts, tolerances, evaluation budgets and bound overrides.


