# what `objective.rs` does
- It is one function that takes a parameter set and returns a single misfit number (the deviance). It runs in these stages.
- Load the observed data on every call. The observed bin edges, frequencies and accuracies come from Flanker_SOA_CSE_ReadData.mat.
- Draw non-decision-time jitter. Each SOA gets 50,000 uniform draws in the range Ter ± std/2.
- Simulate 50,000 trials for each of the 16 conditions. It uses the _pc parameters for sequences starting with a congruent trial (cc, ci) and the _pi parameters for ic and ii. The result is a stored 1500-step path per trial.
- Classify each trial. If the path ever touches −a, the trial counts as an error, even if it touched +a first. Otherwise it is correct if it touches +a, and a timeout if it touches neither. The RT is the step number plus the jitter, and trials with RT ≥ 1500 are discarded.
Compute the predicted accuracy, as correct trials divided by the trials that survive that discarding.
Build the group bin edges by averaging the observed participants' edges, which is what the R script already does. Each simulated correct RT goes into one of the six bins, and anything above the sixth edge goes into none.
Turn observed and predicted counts into "frequencies". Both sides multiply bin counts by an accuracy, and the error bin gets a count times (1 − accuracy).
Scale the predictions to the observed total. Each bin's predicted proportion (over the six correct bins plus the error bin) is multiplied by that condition's observed total frequency. Then the deviance is 2·Σ O·ln(O/E), added up over all bins and conditions. Only the error bins have a guard against zero. The correct bins don't.

