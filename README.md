#  Combining exchangeable p-values

This repository is associated with the paper [Combining exchangeable p-values](https://arxiv.org/abs/2404.03484).

### Abstract

The problem of combining p-values is an old and fundamental one, and the classic assumption of independence is often violated or unverifiable in many applications. There are many well-known rules that can combine a set of arbitrarily dependent p-values (for the same hypothesis) into a single p-value. We show that essentially all these existing rules can be strictly improved when the p-values are exchangeable, or when external randomization is allowed (or both). For example, we derive randomized and/or exchangeable improvements of well known rules like "twice the median" and "twice the average", as well as geometric and harmonic means. Exchangeable p-values are often produced one at a time (for example, under repeated tests involving data splitting), and our rules can combine them sequentially as they are produced, stopping when the combined p-values stabilize. Our work also improves rules for combining arbitrarily dependent p-values, since the latter becomes exchangeable if they are presented to the analyst in a random order. The main technical advance is to show that all existing combination rules can be obtained by calibrating the p-values to e-values (using an α-dependent calibrator), averaging those e-values, converting to a level-α test using Markov's inequality, and finally obtaining p-values by combining this family of tests; the improvements are delivered via recent randomized and exchangeable variants of Markov's inequality.

## Main contents

The script `functions_pvals.R` contains the function introduced in the main paper.

The folder `\Simulations` contains the code to reproduce the simulations in the paper (Section 9 and Appendix F).

