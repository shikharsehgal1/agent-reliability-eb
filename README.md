# Empirical-Bayes reliability testing for stochastic AI-agent tasks

STAT 27855 course project (empirical Bayes, local false discovery rates, multiple
testing). Written in R. Simulation first; one illustrative public dataset.

## Research question

A benchmark has *J* tasks. Task *i* is run *k_i* times and each run succeeds or
fails. With only a few runs per task:

1. Can empirical-Bayes borrowing across tasks improve decisions about which tasks
   meet a prespecified reliability threshold?
2. How misleading is a plug-in local false discovery rate when the mixing
   distribution is not identified by the available run counts?

## Model and notation

For tasks *i = 1, ..., J*:

    theta_i  iid ~ G on [0, 1]
    X_i | theta_i ~ Binomial(k_i, theta_i), independent across tasks
    H_i  : theta_i <= tau   (task is NOT reliable)
    H_i^A: theta_i >  tau   (task is reliable)

    lfdr_i = lfdr(x_i, k_i; G) = int_[0,tau] Bin(x_i; k_i, p) dG(p) / int_[0,1] Bin(x_i; k_i, p) dG(p)

**Two thresholds, two meanings.** `tau` (main value 0.90) is the *reliability
threshold*: the per-run success probability a task must exceed. `q` (main value
0.10) is the *decision threshold*: the cutoff on lfdr (or the BH/BY level) for
declaring a task reliable. Sensitivity: tau in {0.80, 0.90, 0.95}, q in {0.05, 0.10, 0.20}.

A "discovery" means rejecting H_i, i.e. *calling the task reliable*.

The primary analysis uses one binary success criterion per task. Multiple criteria
from the same run are not analysed; they are not independent and treating them as
such would overstate the evidence.

### Assumptions, stated once

* **A1 (exchangeable tasks).** theta_i iid from one G. On a real benchmark the
  tasks are a fixed, curated set; this is a working model, not a sampling fact.
* **A2 (independent runs).** Given theta_i, the k_i runs are independent
  Bernoulli(theta_i). Violated by caching, shared seeds, or correlated infrastructure
  failures (stress case `within_corr`).
* **A3 (independent tasks).** X_i independent across tasks given theta. Violated by
  near-duplicate tasks (stress case `clustered`) or by run-batch effects.
* **A4 (non-informative budgets).** k_i is independent of theta_i. Violated if
  extra runs are spent on tasks that looked borderline.

## The identification problem

For fixed *k*, the marginal of X has k + 1 cell probabilities that depend on G only
through its first k moments. Different G with the same first k moments produce
*exactly the same distribution of observable counts* but can give very different
`Pr(theta <= tau | X = x)`. No amount of data on J (more tasks) fixes this; only
more runs per task do.

The project separates two kinds of uncertainty, which must not be conflated:

| | What it is | What it is not |
|---|---|---|
| **(A) Identified set** | `[inf, sup]` of lfdr(x0, k0; G) over all G (on a fine grid) whose marginal *equals* a given marginal: the true one (population calculation) or the fitted NPMLE marginal (plug-in) | **Not a confidence interval.** It treats the marginal as known. |
| **(B) Inferential interval** | Union of identified sets over a simultaneous 95% confidence box for the marginal (Bonferroni over Clopper-Pearson intervals of the cell frequencies) | Coverage holds if the box covers the true marginal (prob. >= 0.95 under A1-A4), up to grid discretisation. This is the F-localization idea of Ignatiadis & Wager (2022) applied to binomial counts; it is not new here. |

Both are linear-fractional programs over grid weights, solved as LPs after the
Charnes-Cooper transform (`R/bounds.R`).

## Methods compared

| Method | Estimate of Pr(theta_i <= tau \| data) | Decision rule | What the rule means / assumes |
|---|---|---|---|
| Oracle lfdr | true G (simulation only) | lfdr <= q; also step-up | Benchmark. Each call has posterior null prob <= q under the true model. |
| NPMLE EB plug-in | G-hat = Kiefer-Wolfowitz NPMLE on a grid (log-barrier interior point with certified KKT gap; EM/SQUAREM alternative) | lfdr-hat <= q; also step-up | A1-A4 plus "G-hat is close enough to G for this functional", which fails when G is not identified. No frequentist FDR guarantee. |
| Beta-binomial EB | G-hat = Beta(a, b), maximum marginal likelihood, box-constrained, boundary fits flagged | lfdr-hat <= q; also step-up | A1-A4 plus G is Beta. Misspecified for mixtures. |
| Unpooled Bayes | fixed Beta(1, 1) prior for every task; nothing estimated | posterior null prob <= q | A2 plus the prior. No borrowing across tasks. (A prior cannot be fitted from one task's single count.) |
| Exact test + BH | one-sided p_i = Pr(Bin(k_i, tau) >= x_i), least-favourable null theta_i = tau | BH at level q | Valid p-values under A2. BH controls FDR <= q under A3 (independence) or PRDS. |
| Exact test + BY | same p-values | BY at level q | FDR <= q under arbitrary dependence across tasks. Conservative. |
| Box-robust lfdr | U_i = max lfdr(x_i, k_i; G) over all G whose marginal lies in a (1 - q/2) Clopper-Pearson box | U_i <= q/2 (also step-up) | FDR <= q for **any** G under A1-A4 (box fails w.p. <= q/2; otherwise every call has lfdr <= q/2). Averaged over theta ~ G, not conditional on fixed theta. |

"Step-up" = reject the largest set of smallest lfdr values whose mean is <= q, which
targets posterior expected FDP <= q *if the lfdr values are correct*. Neither lfdr
rule controls the frequentist FDR when G-hat is wrong; that is the point of the
simulation.

## Reproducing

Requires R >= 4.3 (developed on R 4.5.2, macOS). From this directory:

```sh
Rscript scripts/install_deps.R   # lpSolve, testthat (+ ggplot2 if missing) into ./.Rlib
sh scripts/fetch_tau_bench.sh    # optional: public tau-bench data, pinned commit
Rscript scripts/run_all.R        # tests, examples, simulation, sensitivity, figures
```

`run_all.R` runs, in order: `tests/run_tests.R`, `scripts/01_examples.R`,
`02_simulation.R` (117 scenarios x 200 replicates; ~30 min on 9 cores; cached per
scenario in `results/raw/`), `03_grid_sensitivity.R`, `04_stability.R`,
`06_tau_bench.R` (if the data were fetched), `05_figures_tables.R`. Set
`NREP=20` for a quick pass. Every replicate has its own seed
(`100000 * scenario_id + replicate`), so results do not depend on parallel scheduling.

## Layout

```
R/            mixing.R       Beta-mixture G, exact oracle lfdr and marginals
              simulate.R     data generation (+ dependence stress cases)
              fit_npmle.R    grid NPMLE: barrier interior point and EM/SQUAREM
              fit_betabinom.R  Beta-binomial EB, unpooled Beta(1,1)
              decisions.R    exact p-values, BH, BY, lfdr threshold and step-up
              bounds.R       identified sets (A) and confidence-box intervals (B) via LP
              metrics.R      FDP, power, near-boundary calls, calibration, Jaccard
              scenarios.R    the factorial design
              replicate.R    one replicate; aggregation with Monte Carlo SEs
scripts/      01..06, run_all.R, install_deps.R, fetch_tau_bench.sh
tests/        testthat suite (formulas, normalisation, NPMLE, LP bounds, reproducibility)
results/      tables/*.csv, figures/*.png, summary.md (all generated)
progress_report.md, references.bib
```

## Simulation design

Factorial, independent model: G (7 shapes) x k (3, 5, 10, 30, heterogeneous
{3,5,10,20}) x J (50, 200, 1000) = 105 scenarios, 200 replicates each, all tau/q
combinations computed on the same replicates. G shapes: single Beta(4,1); separated
and overlapping two-component Beta mixtures with 10% or 40% weight on the reliable
component; mass concentrated near tau; mass near 0 and 1. Pr(theta > 0.9) ranges
0.10-0.56. Stress cases (12 scenarios, reported separately): within-task run
correlation 0.3, and clusters of 10 tasks sharing theta.

k = 30 is included as a contrast: with tau = 0.9, an all-successes task has exact
p-value 0.9^k, which is above 0.10 for every k <= 21, so the exact test cannot make
any call at q = 0.10 in the small-k cells.

## Findings (summary; details and caveats in `progress_report.md`)

* With tau = 0.9, all-successes needs k >= 22 runs before the exact test can reject at
  0.10. Exact test + BH made essentially no reliability calls in the k <= 20 design.
* lfdr is badly non-identified at small k. At k = 3, mixing distributions with
  identical count distributions give lfdr(3|3) anywhere in [0.13, 0.88].
* The NPMLE plug-in makes calls that the oracle would not. Its frequentist FDR
  exceeded q = 0.10 (by more than 2 MCSE) in 10 of 105 primary cells for the
  threshold rule and 31 for the step-up rule, up to 0.20, mostly where G has mass
  near tau. A misspecified Beta-binomial EB reached 0.30, and got worse as J grew.
  When G is well separated from tau, plug-in EB was close to the oracle.
* The plug-in identified set (A) is not a confidence interval: it contained the true
  lfdr in at most 12% of replicates when k >= 10. The confidence-box interval (B)
  covered it in >= 99.5% of replicates, but was nearly [0, 1] unless k = 30.
* A box-robust lfdr rule has provable FDR <= q for any G. Its worst simulated FDR was
  0.004 at q = 0.10, but at tau = 0.9 it made no calls at all unless k = 30.
* On public tau-bench data (k = 4 or 8), no reliability call at tau = 0.9 is
  supported once uncertainty in the marginal is accounted for.

Every table, with Monte Carlo SEs, is in `results/summary.md`.

## Figures (`results/figures/`)

1. `fig1_nonidentification_k3.png`: three G with the same count distribution at k = 3 and different lfdr(3|3).
2. `fig2_identified_sets_by_k.png`: population identified sets for lfdr(k|k) by k and G.
3. `fig3_fdr_by_k.png` / 4. `fig4_power_by_k.png`: FDR and power by k and G (J = 200), +/- 2 MCSE.
5. `fig5_fdr_by_J.png`: FDR vs number of tasks J.
6. `fig6_calibration_k5.png`: reliability diagram of estimated lfdr.
7. `fig7_bounds_coverage.png`: coverage of (A) vs (B).
8. `fig8_selection_curves.png`: Pr(called reliable) as a function of a fixed task's theta.

The cache in `results/raw/` is keyed by scenario id and NREP only, not by code
version. Delete it after changing any code in `R/`.

## Relation to prior work, and what is not new

Nothing here is a new statistical method. Local fdr and empirical-Bayes testing
(Efron et al. 2001; Efron 2010; Sun & Cai 2007), the binomial-mixture NPMLE
(Kiefer & Wolfowitz 1956; Laird 1978; Lindsay 1983, 1995; Wood 1999), the moment
non-identifiability of binomial mixtures (Teicher 1963; Karlin & Studden 1966),
LP bounds on empirical-Bayes functionals and F-localization confidence intervals
(Ignatiadis & Wager 2022), BH/BY (Benjamini & Hochberg 1995; Benjamini & Yekutieli
2001), and pass^k-style repeated-trial agent evaluation (Yao et al. 2024) are all
prior work. The project is an application and a careful simulation study of these
tools in the agent-evaluation setting.

This repository is separate from the author's earlier project `rl-dist-eval`, which
already includes per-task Beta-binomial posteriors with Pr(p > 0.9) reliability
summaries, method-of-moments Beta-binomial EB shrinkage, and BH p-value adjustment.
Those components are not presented as new here. No code was copied from it.

## License

MIT. The tau-bench data (downloaded, not redistributed here) are MIT-licensed by Sierra.
