- rust compiler, fitter and r preprocessing and plotting 
- also: my own decision was not to check against the uploaded .mat files (due to the fact that i had no time, instead i built the error catching in the functions used. 

## step 0: decisions: 
- Scope: CSE models (also CE, if i have time) 
- Corrected: 
    - o.g. code has its quirks: negligible error bin, wrong accuracy reused for ic/ii, error-classification rule:
        - error classification rule: 

## step 1: observed targets: 
- `01_preprocessr.R` -> is output accurate? 
- extend preprocessing to ce (collapsing across the previous trial)

## step 2: rust process skeleton
- on github already
(a beadandoba: csak a vazlat, tenyleg nem kell belemenni a cratekbe stb)

## step 3: the trial simulator: 
- a function that simulates one trial and returns choice and RT
- Don't store whole paths, and stop at the first boundary crossing
- Check it against known results: with constant drift v, boundary a and noise σ = 0.1, accuracy and mean decision time have closed forms. If the simulator doesn't match them, nothing downstream can be trusted

## step 4: the objective (G2 fitting)
- Simulate all conditions in parallel and bin the predicted correct RTs with the observed edges.
- Compute the error mass, rescale to the observed total and calculate the deviance.
- Handle empty bins (the log of zero) explicitly.
- Test it by feeding it data simulated from the model with known parameters. The deviance should be near the noise floor at the true values and grow as you move away.

## step 5: parameter mapping: 
- a table that maps the flat parameter vecor onto each model's full set of parameters (along with the bounds) -> matlab: wrapper files 
- test that every variant gives the right number of free parameters: 
    - CE: 26, 29, 27, 27; CSE: 26 and 30

## step 6: the optimizer: 
- Implement Nelder-Mead with bounds by transformation, the same trick as f`minsearchbnd`. Writing it yourself is about 100 lines. The `argmin` crate is the alternative.
- Make the objective deterministic by fixing the noise seed within a fit.
Add the coarse (8,000 trials) then fine (50,000 trials) stages from the paper, multiple random starts, and a checkpoint file with the best parameters so far.
- Check it first on a simple known function, then on simulated data

## step 7: parameter recovery: 
- Simulate a dataset from chosen parameters, fit it, and see whether you get them back
- With 26 to 30 parameters, some will trade off against each other, and you'll learn which
- This replaces the check against MATLAB's fitted values, which you don't have

## step 8: output for plotting
- Have Rust write the fitted parameters, the deviance, and the predicted quantiles, CDFs and accuracies per bin to CSV

## step 9: plots (R)
- CDFs, delta plots, CAFs, observed against predicted

## step 10: model comparison
- small r script that takes the deviances and parameter counts and runs tha chi-squared difference tests 

## step 11: real fits: 
- run full pipeline 


# error handling, changes in my code compared to the original: 
- debug functions: 
    - rust reads only the csv files that the r code writes. 
- error bin: proper 7th bin WITH the sim error proportion x the observed total
- ic/ii accuracy: each condition uses its own simulated accuracy
## checks that don't need the binary `.mat` files
- closed form sim checker (CAF, CDF, delta plots)
- -> check that is the most convincing: a fit that reproduces the observed pattern
- parameter recovery
- check parameter estimates reported in the paper

