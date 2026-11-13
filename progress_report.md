# Progress report: empirical-Bayes reliability testing for AI-agent tasks

STAT 27855 project. Status as of 2026-11-12. All numbers below come from
`Rscript scripts/run_all.R` as committed (tests passing; 117 scenarios x 200
replicates); full tables with Monte Carlo SEs are in `results/summary.md`.

## 1. Question and setup

Tasks i = 1..J have per-run success probability theta_i ~ G (iid), and
X_i | theta_i ~ Bin(k_i, theta_i). A task is *reliable* if theta_i > tau. We
call a task reliable when its posterior null probability lfdr_i = Pr(theta_i <=
tau | X_i) is at most q, where tau = 0.90 is the reliability threshold and
q = 0.10 the decision threshold. Questions: (i) does empirical-Bayes borrowing across
tasks improve these calls when k is small, and (ii) how misleading is the plug-in
lfdr when G is not identified? With k runs, the count distribution depends on G only
through its first k moments, so lfdr is generally not identified.

## 2. What has been done

* **R analysis repository** (`R/`, `scripts/`, `tests/`), separate from the
  author's earlier `rl-dist-eval` project, which already contained Beta-binomial
  posteriors, MoM empirical-Bayes shrinkage and BH adjustment; those are not claimed
  as new here.
* **Seven methods**: box-robust lfdr rule with guaranteed FDR (Section 3); oracle lfdr (true G); grid NPMLE plug-in (Kiefer-Wolfowitz;
  log-barrier interior point with a certified KKT gap, plus EM/SQUAREM); Beta-binomial
  EB (bounded MLE with boundary flags); unpooled Beta(1,1); exact one-sided binomial
  test at theta = tau with BH and BY. lfdr rules are run as a threshold rule
  (lfdr <= q) and as a step-up rule (mean lfdr of the called set <= q).
* **Identification bounds** by linear programming (Charnes-Cooper transform):
  (A) the identified set of lfdr(x0, k0) under a *given* marginal, which is not a
  confidence interval; and (B) the union of identified sets over a simultaneous 95%
  Bonferroni Clopper-Pearson box for the marginal. (B) is an inferential interval,
  an instance of the F-localization approach of Ignatiadis & Wager (2022).
* **Simulation**: 7 shapes of G x k in {3, 5, 10, 30, heterogeneous 3/5/10/20} x J in
  {50, 200, 1000}, 200 replicates each, all tau in {0.8, 0.9, 0.95} and q in {0.05,
  0.1, 0.2} computed on the same replicates; 12 dependence stress cells. Also: NPMLE
  grid/solver/smoothing sensitivity (20 datasets x 3 cases x 9 configurations) and
  discovery-set stability with tasks held fixed (100 redraws).
* **Public data**: tau-bench historical trajectories (Yao et al. 2024; MIT
  license; pinned commit 59a200c). They contain task_id, trial and reward in {0, 1}:
  airline (50 tasks) and retail (115 tasks), with 4 trials for GPT-4o and 8 for
  Claude 3.5 Sonnet (new). This is a secondary illustration only (Section 4).
* **Tests**: 24 testthat tests, 140 expectations, including checks of the robust
  rule's monotone early stopping and of its FDR in a hard case. They check the closed-form lfdr
  against numerical integration, normalisation, BH/BY against `p.adjust`, the NPMLE
  optimality certificate, marginal recovery, the LP bounds against an analytic k = 1
  example, nesting of (A) inside (B), detection of an infeasible marginal, and seed
  reproducibility. All pass.

## 3. Main results (tau = 0.90, q = 0.10 unless stated)

**All successes is weak evidence.** The exact p-value for k/k is 0.9^k, so the test
cannot reject at 0.10 until k = 22 (at 0.05, k = 29). Under a Beta(1,1) prior,
3/3 gives Pr(theta > 0.9) = 0.34 and 10/10 gives 0.69. Consequently **exact test + BH
made essentially no reliability calls anywhere in the tau = 0.9 design**. Its power
was 0 except with k = 30 and G having mass near 0 and 1, where it was at most 0.19. At tau = 0.8, k = 30 it worked as advertised: power 0.13-0.89 with
FDR <= 0.005. It is conservative, because the p-values are discrete and evaluated at
the least-favourable null.

**Non-identification is large at small k (Fig. 1-2).** For the overlapping mixture
at k = 3, LP finds mixing distributions with *the same* count distribution
(to 1e-9) giving lfdr(3|3) = 0.13 and 0.88; the true value is 0.46. Pr(theta > 0.9)
itself ranges over [0.08, 0.75]. The identified set for lfdr(k|k) shrinks with k but
stays wide when G has mass near tau: [0.01, 0.92] at k = 10 and [0.04, 0.68] at
k = 20 for the near-tau G. **Caveat:** the extremes are attained by placing an atom
exactly at tau versus just above it, so part of the width reflects the
discontinuity of Pr(theta <= tau) at a sharp threshold, not only moment
indeterminacy.

**Plug-in EB does not control frequentist FDR, and fails where it matters (Fig. 3-5).**
Across the 105 primary cells (J and k varied), the NPMLE threshold rule had FDR
significantly above 0.10 (mean - 2 MCSE > 0.10) in 10 cells and the NPMLE step-up rule
in 31. The worst cell was G with mass near tau and heterogeneous k: FDR 0.20 (MCSE 0.02),
while the oracle made no calls at all. The Beta-binomial EB, misspecified for that G,
reached FDR 0.30 at J = 1000, and its error *grew* with J (0.20 at J = 50), as
expected for a bias that more tasks cannot remove. Even the fixed Beta(1,1)
prior reached FDR 0.23 (near-tau G, k = 30). When G is well separated from tau,
the plug-in did well: e.g. separated G with 40% reliable component, k = 10, J = 200:
NPMLE FDR 0.018 and power 0.84, versus oracle 0.008 and 0.79.

**Calibration (Fig. 6).** The NPMLE lfdr is overconfident exactly in the decision
region. With k = 5, J = 200 and overlapping or near-tau G, tasks with estimated lfdr
in [0, 0.3] were null 39-69% of the time. The mean |lfdr-hat - lfdr| at J = 200
ranged from 0.004 (separated, k = 30) to 0.27 (near-tau, k = 10).

**Instability (Fig. 8; grid study).** With tasks fixed and only runs redrawn
(overlapping G, k = 10), the NPMLE called 27 tasks reliable on average, 8 of them
falsely; the oracle called none. The mean Jaccard overlap between discovery sets was
0.46. In the grid study (near-tau G, k = 10, J = 200), the NPMLE lfdr(10|10) had SD
0.38 across datasets around a true value of 0.39. Grid
size (51 to 801 points) and solver/start (barrier vs EM) changed lfdr by about 0.01-0.06;
Dirichlet smoothing changed it by 0.2-0.5 and removed most calls. At equal
log-likelihood the fitted G is not unique, and different solvers return different
maximisers.

**Identified set vs inferential interval (Fig. 7).** For lfdr at all-successes,
the plug-in identified set (A) under the fitted marginal often *excluded* the true
value: coverage at most 12% in every cell with k >= 10. The fitted NPMLE marginal sits on
the boundary of the moment space, so (A) collapses to a near-point. This is the
empirical reason (A) must not be read as a confidence interval. The confidence-box
interval (B) contained the oracle lfdr in at least 99.5% of replicates in every cell
(nominal 95%; conservative). Its mean width exceeded 0.9 for k <= 10 except when G had little
mass near tau. Even then (edges or separated G, J >= 200) it narrowed only to
0.17-0.85 and never certified lfdr <= q. It certified lfdr(k|k) <= 0.10 only with
k = 30 and edges or separated G: in 81-100% of replicates for edges G (J >= 200),
separated G with 40% reliable component (J >= 200) and with 10% (J = 1000); in 7.5%
for edges G at J = 50; never otherwise.

**A rule with a guarantee: box-robust lfdr.** Build the (B) confidence box at level
1 - q/2, let U_i be the largest lfdr(x_i, k_i; G) over all G consistent with it, and
call task i reliable when U_i <= q/2. On the event that the box covers the true
marginal, every call has true lfdr <= q/2, so E[V/R | X] <= q/2. The box fails with
probability <= q/2. Hence **FDR <= q for any G** under A1-A4, without assuming a
shape for G or consistency of an estimator. Because lfdr(x, k; G) is nonincreasing
in x for every G, only a few LPs per run budget are needed. In simulation the
guarantee held with a wide margin: maximum FDR 0.004 over all 105 primary cells at
q = 0.10, against 0.20 for the NPMLE plug-in and 0.30 for the Beta-binomial. The
cost is power. At tau = 0.9 it made calls only with k = 30 and G having little mass
near tau: power 0.40-0.69 (edges G, J >= 200) and up to 0.55 (separated G, J = 1000).
With k <= 10 it made **no calls in any cell**. At tau = 0.8 it sometimes had power where
exact test + BH had none, through borrowing across tasks: edges G, k = 10, J = 200,
power 0.23 vs 0; heterogeneous k, 0.18 vs 0.06. With k = 30, BH was usually more
powerful. The two guarantees differ: BH's holds for fixed theta's, the box-robust
rule's is averaged over theta ~ G.

**Dependence stress.** Within-task run correlation 0.3 (runs not independent) raised
FDR in the overlapping-G cells to 0.26-0.31 (NPMLE) and 0.42 (Beta-binomial, k = 10).
Clusters of 10 tasks sharing theta changed FDR little at J = 200.

## 4. Public-data illustration (tau-bench)

On Claude 3.5 Sonnet (new) retail (J = 115, k = 8; 41 tasks at 8/8), the NPMLE
plug-in gives lfdr(8|8) = 0.048 and so calls 41 tasks reliable at tau = 0.9. The
Beta-binomial EB gives 0.13 and calls none, as do the unpooled and BH rules. For
GPT-4o retail (k = 4), the plug-in lfdr(4|4) is 0.0008 at tau = 0.8 but 0.9996 at
tau = 0.9. The 95% interval (B) for lfdr(k|k) is about [0, 1] in every file and
threshold except one ([0.00, 0.47]). So these data cannot support any reliability
call at tau = 0.9, which matches the simulation. A chi-square test found no
trial-batch effect (p = 0.39-0.98), but conditional independence of trials and the
exchangeability of a curated task set remain unverified assumptions. The decoding
settings used to generate these files are not documented.

## 5. Difficulties and how they were handled

* EM with SQUAREM hit its iteration cap without a small KKT gap. I replaced it as
  the default with a log-barrier Newton method whose gap is certified. 98.6% of fits
  met the 1e-6 gap flag, and every fit had gap <= 1.07e-6.
* lpSolve occasionally fails numerically (status 3 or 5) on badly scaled problems,
  e.g. k = 50 with lfdr near 1e-5. The code retries with other scaling modes;
  remaining failures (at most 4 of 400 LPs in a cell) are reported as NA and
  excluded from the coverage denominators.
* An early test assumed NPMLE lfdr recovery at k = 60. It failed (errors up to 0.3
  near tau, with sign changing across seeds), so the test was changed to what the
  theory supports: marginal recovery, and lfdr recovery only for G away from tau.
* A filtering bug that double-counted the independent control cells of the stress
  design in the primary tables was found and fixed before these results were produced.

## 6. Limitations and remaining work

* (B) uses a Bonferroni box and is very conservative. A KS- or
  likelihood-ratio-based localisation, as in Ignatiadis & Wager, could be tighter.
  Coverage was checked only for x0 in {k0, k0 - 1}, and is approximate because G is
  restricted to a grid.
* The box-robust rule controls FDR but is far below its nominal level (0.004 vs
  0.10), and has no power at k <= 10. Next: tighten the box (KS or
  likelihood-ratio localisation, a better split of q between box and threshold), and
  test the rule under the dependence stress cases at k = 30, where it actually makes
  calls (current stress cells use k <= 10, so its zero FDR there is vacuous).
* A margin-of-indifference version of the hypothesis (e.g. theta > tau + delta)
  would remove the boundary-atom artefact in the bounds; not yet done.
* Real data: tau-bench is small (J <= 115, k <= 8) and is used only as an
  illustration; no claim is made about the models evaluated.
* Literature: the methods are standard (Efron; Kiefer-Wolfowitz; Ignatiadis &
  Wager). I have not done a systematic review establishing a gap, and nothing here is
  claimed as a methodological contribution.
