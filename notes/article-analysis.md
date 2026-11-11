## exclusion criteria: 
- Trials with RTs faster than 200 ms and slower than 1500 ms were discarded from data analysis
- first trial of each block 
- error trials, and trials following error trials
-> discarding around 10% of all trials -> no exact number ??

## analysis: 
- dv: mean RT 
- rmANOVAs: three within-subject factors: 2 (prev-cong) x 2 (current-cong) x 4 (SOA: -100, -50, 0, +100)
- mean error rate did not exceed 6% in any SOA condition

### accuracy rates: table from Appendix A. (Table 6, 7)
a. Table 6. Stroop task (in terms of the three-way interaction)
SOA     | cC    | cI    | iC    | iI 
------- | ----- | ----- | ----- | -----
100-ms  | 1.00  | .97   | 1.00  | .98
0-ms    | .99   | .96   | .99   | .98
-50-ms  | 1.00  | .96   | .99   | .98
-100-ms | .99   | .94   | .99   | .97

b. Table 7. Flanker task (in terms of the tree-way interaction)
SOA     | cC    | cI    | iC    | iI 
------- | ----- | ----- | ----- | -----
100-ms  | 1.00  | .98   | 1.00  | .98
0-ms    | 1.00  | .97   | 1.00  | .97
-50-ms  | .99   | .97   | 1.00  | .98
-100-ms | 1.00  | .94   | 1.00  | .97

## Full ANOVA analysis results for the Flanker task: 
Factor                         | df    | F      | p     | partial eta squared 
------------------------------ | ----- | ------ | ----- | --------------------
curr.congruency                | 1,39  | 542.85 | >.001 | .93
previous congruency            | 2,39  | 0.01   | .943  | <.01
SOAs                           | 3,117 | 5.27   | .014  | .12
curr.cong x prev.cong.         | 1,39  | 2.60   | .115  | .84
curr.cong x SOAs               | 3,117 | 117.16 | >.001 | .75
prev.cong. x SOAs              | 3,117 | .70    | .556  | .02
prev.cong x curr.cong. x SOAs  | 3,117 | 6.21   | .001  | .14

- main effect of *currect congruency* 
    - RT on incong trials = 526 ms
    - RT on cong trials = 464 ms
    - CE = 526-464 = 62 ms
- main effect of *SOAs*
    - mean RT was significantly faster on sequential presentation: asynchronous onsets of target and distractor; 495, 494 ms, 483 ms in the order of SOA = -100, -50, and 100 ms v.s 505 ms in the simultaneous presentation
- two-way interactions: *curr.cong x SOAs* [Figure 6. - CE vs SOA plot, color = empirical vs. prediction]
    - Congruency effects increased in negative SOAs compared to the simultaneous condition 
    - -100 ms: 89 ms
    - -50 ms: 73 ms
    - 63 ms in the 0-ms SOA
- *prev.cong x curr.cong. x SOAs*
    - TEST OF SIMPLE EFFECTS: significant small CSEs only:
        - in the -50ms condition: F(1,39) = 4.79, p = .035, partial eta sq. = .11
        - in the -100ms condition: F(1,39) = 20.57, p < .001, partial eta sq. = .35
    - Fig. 7., 8. other conditions: numerically reversed or no CSEs (0ms, 100ms)

## Full ANOVA analysis results for the Spatial Stroop task:
Factor                         | df    | F      | p     | partial eta squared 
------------------------------ | ----- | ------ | ----- | --------------------
curr.congruency                | 1,39  | 314.39 | <.001 | .89
previous congruency            | 1,39  | 24.46  | <.001 | .39
SOAs                           | 3,117 | 37.16  | <.001 | .49
curr.cong x prev.cong.         | 1,39  | 95.46  | <.001 | .71
curr.cong x SOAs               | 3,117 | 12.09  | <.001 | .24
prev.cong. x SOAs              | 3,117 | .80    | .496  | .02
prev.cong x curr.cong. x SOAs  | 3,117 | 4.43   | .005  | .10

- *current cong*
    - mean RT on incong trials = 518 ms
    - mean RT on cong trials = 463 ms
    - CE = 54 ms
- *prev-cong*
    - mean RT on prev-incong trials = 488 ms
    - mean RT on prev-cong trials = 492 ms
- *SOA*
    - mean RT was generally faster in longer SOA conditions
    - -100ms SOA => 474 ms
    - -50ms SOA => 480 ms
    - 0ms SOA => 506 ms
    - 100ms SOA => 499 ms
- *prev-cong x current cong*
    - larger congruency effect in trials following a congruent trial than following an incongruent trial (CSE = 25 ms)
    - TEST OF SIMPLE EFFECTS: 
        - mean RT was faster in congruent trials following a congruent trial (cC; 454 ms) vs. congruent trials following an incongruent trial (iC; 472 ms) trials (F(1,39) = 104.76, p<.001, partial eta sq. = .73), and faster in incongruent trials following another incongruent trial (iI; 513 ms) than incongruent trials following a congruent trial (cI; 522 ms) trials (F(1,39) = 38.40, p<.001, partial eta sq. = .50).
- *current cong x SOA*
    - smaller CE as the distractor leads the target for longer durations
    - smaller CE for the -50ms and the -100ms SOAs compared to the 0ms SOA and a larger CE for the 100ms SOA compared to the 0ms SOA
    [Fig. 3.]
- *prev-cong x SOA*
    - post-conflict slowing being larger at -50-ms SOA (6 ms) and weakest when the target leads (+100-ms SOA; 2 ms)
- *prev-cong x current cong x SOA*
    - the CSE increases as the distractor leads for a longer time
    - -100ms CSE = 34 ms
    - -50ms CSE = 22 ms 
    - 0ms CSE = 23 ms
    - +100ms CSE = 20 ms




#### Altogether pattern regarding the results: 
The behavioral data showed two main patterns. Firstly, as the distractor lead-time increased, the congruency effect decreased in the spatial Stroop task while it increased in the flanker task. Secondly, the CSE increased in the spatial Stroop task with a larger temporal lead-time for the distractor, yet it was mostly absent across SOAs in the flanker task. In the next section, we will use RDMC to investigate the mechanisms underpinning the variations of congruency effect and CSEs across SOAs.


