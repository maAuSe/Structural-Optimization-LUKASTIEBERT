# Questions for Prof. Schevenels — meeting on 11 May 2026

## Context

The assignment guidelines state that every group is expected to request guidance at least once *before* including the results in the report. The report deadline is 15 May 2026 and the presentation is on 22 May 2026, so this meeting is our last opportunity to align on the methodological choices that shape the headline conclusions of the report. The two questions below are therefore the ones we cannot easily revisit after submission. Smaller follow-ups are listed at the end.

For reference, the current state of our work is:

- All four assignment steps have been implemented and produce reproducible results (see `matlab/run_assignment.m`, `matlab/results/*.mat`, and the figures in `matlab/figures/`).
- Sensitivity verification against central finite differences passes with relative errors of $5.7\times 10^{-7}$ (classical SIMP), $1.1\times 10^{-5}$ (MTOP per-cell), and $5.4\times 10^{-6}$ (full filter–Heaviside–MTOP chain).
- The report draft (`report/main.pdf`) is structured around a $9$-run experiment matrix; numbers cited below come directly from the result summary files in `matlab/results/`.

---

## Question 1 — Convergence handling for the density-filtered runs

### What we observe

In every run that uses density filtering — classical OC (step 3), MTOP OC (step 3), and MTOP MMA (step 4) — the standard design-change tolerance $\max_i |\Delta x_i| < \varepsilon = 0.01$ is **never satisfied**. Instead the maximum absolute design change plateaus between roughly $0.04$ and $0.07$ for the rest of the run, while the compliance has long since stabilized to within $0.1\,\%$ of its final value. Concretely:

| Run                               | Iter. | Final $\max|\Delta x|$ | Stop reason     |
|-----------------------------------|------:|-----------------------:|-----------------|
| Classical OC, density filter      |   237 |                  ~0.05 | compl. plateau  |
| MTOP OC, density filter           |   299 |                 0.0516 | compl. plateau  |
| MTOP MMA, density filter          |   419 |                  ~0.05 | compl. plateau  |

To stop these runs cleanly we added a second termination criterion alongside the design-change one: the relative spread of the compliance over a sliding window of 20 iterations must drop below $10^{-4}$. With this criterion in place all three runs terminate with a well-defined compliance and a residual high-frequency oscillation on the boundary of the structural members.

We interpret the oscillation as the well-known incompatibility between the multiplicative OC update and density filtering on a fine mesh: the filter spreads each design-variable change over a neighbourhood of cells, but the OC update rescales each variable individually with no built-in damping, so high-frequency components of the design alternate sign from one iteration to the next even though the compliance has already converged. The same story plays out under MMA — switching the optimizer alone does *not* fix the change tolerance — and only the Heaviside projection ($\eta\in\{0.3, 0.5, 0.7\}$) produces designs that genuinely satisfy $\max|\Delta x|<\varepsilon$ once $\beta$ has reached $\beta_{\max}=16$.

### What we want to ask

Is the compliance-plateau termination criterion the framing you want for the report, or would you prefer that we resolve the oscillation directly with one of the textbook remedies before we hand in? The candidates we have considered are:

1. **Smaller OC move limit** (e.g. $m=0.05$ instead of $m=0.2$): would damp the oscillation at the cost of more iterations and therefore a larger wall-clock budget for the density-filtered cases.
2. **Penalty continuation on $p$** (e.g. ramp from $p=1$ to $p=3$ over the first 50–100 iterations): standard SIMP technique that is also expected to mitigate post-plateau oscillation.
3. **Explicit damping in the OC update** (e.g. $x^{\mathrm{new}}_i \leftarrow (1-\alpha)\,x_i + \alpha\, x^{\mathrm{OC}}_i$ with $\alpha < 1$): direct fix for the multiplicative-update issue, but not part of the lecture form of OC.
4. **Keep the compliance-plateau criterion**, document it, and report it as the genuine convergence behaviour of OC + density filter on a fine mesh — i.e. treat the oscillation as a result, not a problem to be fixed.

We would prefer option 4 because we believe the oscillation *is* the result here (and the comparison between formulations stays clean), but we are conscious that introducing a non-standard convergence criterion in 3 of 9 runs may look ad hoc. Your guidance would help us frame Section 7.4 ("Density filtering with vanilla optimality criteria") of the report.

---

## Question 2 — Reporting the cost/accuracy trade-off of MTOP

### What we observe

Every MTOP design has two compliance values:

- the **native compliance**, computed inside the MTOP formulation on the coarse $120\times 40$ analysis mesh with the multiresolution stiffness matrix; and
- the **fine-mesh compliance**, obtained by re-assembling the same density field on the $600\times 200$ classical FE model and solving once.

The two values differ by roughly $2.5\,\%$ for the sensitivity-filtered case ($219.069$ vs. $224.696$) and by about $2.3\,\%$ for the density-filtered case ($235.363$ vs. $240.984$). This gap is *not* an optimization error. It is the residual quadrature error of the per-cell midpoint integration in Eq. (5) of the report relative to the analytical Q4 stiffness, and it scales as $\mathcal{O}(h^2)$ with $h = 1/\sqrt{N_{\mathrm{c}}}$.

Our current convention in the report is:

- **Layout comparison** (classical vs. MTOP): we use the fine-mesh compliance for the MTOP design. This gives a gap of $0.045\,\%$ for sensitivity filtering and $0.04\,\%$ for density filtering, well below the SIMP modeling error itself. We argue this is the only fair way to compare the *mechanical* quality of the two designs.
- **Speedup comparison**: we report a wall-clock speedup of **$21.4\times$** for sensitivity filtering and **$5.07\times$** for density filtering. The fine-mesh re-analysis time is *not* charged to the MTOP run; we treat it as a one-off verification step rather than as part of the optimization budget.
- **Headline interpretation** (Section 7.3 "Origin of the speedup"): the per-phase wall-clock breakdown shows the FE solve is $40\times$ faster on the coarse mesh, but the overall speedup is only $21.4\times$ because the cone-kernel filter and the OC bisection scale with the density mesh and are mesh-independent. The speedup gap between sensitivity ($21.4\times$) and density filtering ($5.07\times$) is fully explained by extra cone convolutions in the density-filter chain rule and inside the OC bisection.

### What we want to ask

We have three sub-questions, all about the same underlying issue:

**2(a) — Which compliance is "the result" of MTOP?** Should we present the fine-mesh value as *the* compliance of the MTOP design (treating the native value as a verification check), or do you consider the native MTOP compliance the legitimate result and the gap a known feature of the formulation? Our current draft uses the fine-mesh value because it allows direct comparison with the classical reference, but we are aware that a user who only ever runs MTOP would never pay the cost of the fine-mesh re-analysis and would simply report the native value.

**2(b) — Should the fine-mesh re-analysis cost be charged to MTOP in the speedup figure?** If the answer to 2(a) is "fine-mesh", consistency would suggest including the re-analysis time in the MTOP wall-clock budget. For the sensitivity-filtered case this would add roughly the cost of one classical FE solve (~3 s on our machine) to the $12.28$ s of MTOP, dropping the speedup from $21.4\times$ to about $17\times$. We do not currently do this, and we want to confirm this is the right framing.

**2(c) — Is the headline speedup the right narrative, or should we extend the experiment matrix?** The Nguyen et al. paper does not report timings, so the speedup is the genuine novelty of our assignment. The per-phase breakdown already explains *why* the speedup is what it is, but we could strengthen the conclusion with one extra parameter sweep, for example:

- different number of density cells per element ($N_c \in \{3\times 3, 5\times 5, 7\times 7\}$), to show how the speedup vs. accuracy trade-off depends on $N_c$;
- a larger mesh (e.g. $1200 \times 400$ density cells), to show whether the speedup grows with problem size as the asymptotic Cholesky scaling argument would predict;
- a different filter radius, to show the influence on the per-density-cell phase.

We have time to fit in *one* such sweep before the deadline. Which (if any) would you find most informative?

---

## Smaller follow-ups (only if there is time)

1. **Group composition.** The assignment specifies "groups of three students" but our group is two (Lukas Campaert and Tiebert Lefebure). We have been working as a pair from the start; we want to confirm there is no issue with submitting as a group of two.

2. **GenAI transparency.** Section 9 of our report describes how we used GenAI in three roles (repository scaffolding, translation of the lecture examples into production runners, drafting of the multiresolution stiffness assembly and the verification routine). Is this level of granularity what you expect, or would you prefer a more itemized statement (per chapter, per figure, etc.)?

3. **Presentation focus.** The presentation is 10 minutes plus 10 minutes of questions. Given the audience (other students who took the course), we plan to spend most of the time on the per-phase wall-clock breakdown and the Heaviside results, and treat the convergence issue of question 1 as a single slide. Does that match what you would expect?
