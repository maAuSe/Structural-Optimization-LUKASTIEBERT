# Questions for Prof. Schevenels


## Question 1 - Convergence handling for the density-filtered runs

### What we observe

In our current density-filtered production runs -- classical OC (step 3), MTOP OC (step 3), and MTOP MMA (step 4) -- the standard design-change tolerance
$$
\max_i |\Delta x_i| < \varepsilon = 0.01
$$
does not become the active stopping condition. The compliance stabilizes much earlier, while the maximum absolute design-variable change remains at the level of a small residual oscillation.

| Run                          | Iter. | Final $\max|\Delta x|$ | Stop reason        |
|------------------------------|------:|-----------------------:|--------------------|
| Classical OC, density filter |   237 |                 0.0491 | compliance plateau |
| MTOP OC, density filter      |   299 |                 0.0516 | compliance plateau |
| MTOP MMA, density filter     |   419 |                 0.0470 | compliance plateau |

To stop these runs consistently, we added a second termination criterion alongside the design-change one: the relative spread of the compliance over a sliding window of 20 iterations must drop below $10^{-4}$. With this criterion, all three density-filtered runs terminate with a stable compliance value and a residual high-frequency change concentrated near the member boundaries.

Our interpretation is deliberately cautious. For the OC runs, the behavior appears consistent with the known difficulty of combining density filtering with the multiplicative OC update: the filter spreads each design-variable change over a neighborhood of cells, while the OC update still rescales variables individually. That can leave high-frequency components in the design variables even after the physical response has stabilized. The MMA density-filtered run shows a similar residual-change issue under our current settings, so switching optimizer alone did not remove the phenomenon, but we do not want to overclaim that the mechanism is identical.

### What we want to ask

Would you consider the compliance-plateau criterion an acceptable way to frame these density-filtered results, or should we resolve the residual oscillation before submission?

The options we see are:

1. **Keep the compliance-plateau criterion**, document it explicitly, and report the residual design-variable oscillation as part of the observed behavior of the density-filtered formulation.
2. **Use a different convergence measure**, for example the change in physical density $\rho=\mathcal{F}(x)$ rather than the raw design-variable change. This may better reflect the quantity entering the FE analysis, but it would differ from the lecture stopping criterion.
3. **Apply an OC-specific remedy**, such as a smaller OC move limit ($m=0.05$ instead of $m=0.2$) or explicit damping in the OC update. This would likely damp the OC oscillation, but it changes the lecture form of OC and does not directly address the MMA density-filtered run.
4. **Add continuation or stronger regularization**, for example introducing a new SIMP-penalty continuation that starts at $p=1$ and ramps to the current assignment value $p=3$. This is more standard than ad hoc damping, but it changes the experiment matrix and would require rerunning the affected cases.

Our current preference is option 1, because it keeps the comparison between formulations clean and treats the residual oscillation as a convergence feature that should be reported rather than hidden. We would like to know whether you agree, or whether you would expect one of the remedies above before the report is submitted.

---

## Question 2 - Reporting the cost/accuracy trade-off of MTOP

### What we observe

Every MTOP design has two compliance values:

- the **native compliance**, computed inside the MTOP formulation on the coarse $120\times 40$ analysis mesh with the multiresolution stiffness matrix; and
- the **fine-mesh compliance**, obtained by re-assembling the same final density field on the $600\times 200$ classical FE model and solving once.

The two values differ by roughly $2.5\,\%$ for the sensitivity-filtered case ($219.069$ vs. $224.696$) and by about $2.3\,\%$ for the density-filtered case ($235.363$ vs. $240.984$). We do not interpret this as an optimization failure: the layouts agree closely with the classical reference, and the sensitivities verify against finite differences. Our current explanation is that the gap is the residual approximation/integration error of the per-cell midpoint integration in Eq. (5) of the report relative to the analytical Q4 stiffness. For smooth integrands, this error is expected to scale as $\mathcal{O}(h^2)$ with $h=1/\sqrt{N_{\mathrm{c}}}$, but we have not yet run an $N_{\mathrm{c}}$ sweep to demonstrate that scaling directly.

Our current convention in the report is:

- **Layout/mechanical comparison** (classical vs. MTOP): we use the fine-mesh compliance for the MTOP design. This gives a gap of $0.045\,\%$ for sensitivity filtering and about $0.04\,\%$ for density filtering, so the mechanical quality of the final layouts is compared on the same FE model.
- **Optimization-cost comparison**: we report wall-clock speedups of **$21.4\times$** for sensitivity filtering and **$5.07\times$** for density filtering. The fine-mesh re-analysis time is not charged to the MTOP optimization run; we treat it as a one-off verification solve.
- **Speedup narrative**: the per-phase wall-clock breakdown shows that the FE solve is about $40\times$ faster on the coarse analysis mesh, but the overall sensitivity-filter speedup is only $21.4\times$ because filtering, multiresolution assembly, and the OC update still include operations that scale with the density mesh. The smaller density-filter speedup ($5.07\times$) is largely explained by the extra cone-kernel operations in the density-filter chain rule and volume checks.

### What we want to ask

We have three related questions about how to present this trade-off.

**2(a) - Which compliance should be the primary MTOP result?** Should the report present the fine-mesh compliance as the main compliance of the final MTOP design, with the native MTOP value reported as the optimization objective used during the run? Or would you prefer the native MTOP compliance to remain the main result, with the fine-mesh value treated only as a cross-check against the classical formulation?

Our current draft uses both values, but uses the fine-mesh value for direct layout-quality comparisons because it puts classical and MTOP designs on the same analysis model.

**2(b) - Should the fine-mesh re-analysis cost be charged to MTOP?** If the fine-mesh compliance is used as the main comparison value, should the cost of that final verification solve be included in the MTOP wall-clock budget? For the sensitivity-filtered case, this would add roughly one classical FE solve to the $12.28$ s MTOP runtime and reduce the headline speedup from $21.4\times$ to about $17\times$.

Our current framing separates the two costs: MTOP optimization time for speedup, and fine-mesh re-analysis as an optional verification step. We would like to confirm whether that distinction is acceptable.


---

## Smaller follow-ups (only if there is time)

3. **Presentation focus.** The presentation is 10 minutes plus 10 minutes of questions. Given the audience (other students who took the course), we plan to focus on the per-phase wall-clock breakdown and the Heaviside results, and treat the convergence issue from question 1 as one short slide. Does that match what you would expect?


### FEEDBACK Prof. Schevenels

1) Sensitivity verification (Section 5).
Add plots & figures. Look at additional PDF with finite-difference approach. Compute sensitivities for a couple of design variables, put them in a comprehensive table. Choose the best method of the 3 presented in the PDF. Or use all methods and compare results.
FILE: "Sensitivity-verification-with-finite-differences.pdf"

2) Experiment matrix (Section 6.1).
Change the name to e.g. "Experimental setup".

3) Density filtering results (Section 7.3 & 7.4).
Compliance plateau convergence criterion instead of classical stop criterion. Add figures of sensitivities for additional proof.

4) Runtime comparison (Section 7.7).
Sort cases differently. Put "sensitivity filtering", "density filtering" and "Heaviside projection" together (per 3). Sensitivity filtering and density filtering are actually different problems.
Native compliance vs. fine-mesh compliance comparison. Fine-mesh compliance should be used as the main comparison value.

5) Coarse-coarse compliance.
New analysis mesh (coarse). Compare compliance values. Is MTOP really that necessary? Coarser mesh with a grayer design is often sufficient.