# MTOP oral defense: anticipated questions and answers

Preparation document for the oral defense of the *Structural Optimization* assignment (B-KUL-H0T88A), on multiresolution topology optimization (MTOP) of the half MBB beam.

Twenty questions a technically-minded examiner is likely to ask, with model answers grounded in the report and the course material. Both presenters should be able to answer every question. The answers are written in full for study; when spoken, compress each to roughly 30 to 60 seconds and lead with the one-sentence core, then justify.

## A. The MTOP method

### Q1. What exactly does MTOP approximate, and why is a fine-resolution design still valid when the analysis runs on a coarse mesh?

MTOP separates two meshes that are normally identical. The design (density) field stays at 600 by 200; the displacement field is discretized on the coarse 120 by 40 mesh. So the only thing approximated is the displacements: they are piecewise-bilinear over each coarse element, meaning the analysis cannot resolve displacement variation finer than one coarse element. The material distribution is still resolved at the full 600 by 200.

Inside a coarse element the 25 density cells enter only through the element stiffness K_e = Σ_k [E_min + ρ_ek^p (E0 − E_min)] I_k. They modulate stiffness pointwise, but the strain field they multiply is the smooth coarse-element strain field. This is valid as long as the true displacement solution is itself smooth over a coarse element, which holds when the minimum feature size (set by the filter, radius 24 density cells, about 5 coarse elements) is comfortably larger than one coarse element. We then re-analyze the final design on the full 600 by 200 mesh to confirm the coarse displacement assumption did not distort the compliance: it differs by only 0.045 %.

### Q2. Explain the midpoint-quadrature element stiffness. Why 25 sub-cell templates, and how accurate is it versus exact Q4 integration?

A standard Q4 element stiffness is K0 = integral of B^T D0 B over the element, normally evaluated with 2 by 2 Gauss quadrature. In MTOP the material factor E(ρ) varies inside the element (one value per density cell), so we integrate with one quadrature point at the centroid of each of the 25 density cells: K_e = Σ_k [E_min + ρ_ek^p (E0 − E_min)] I_k, with I_k = A_k B^T(ξ_k) D0 B(ξ_k) and A_k = 1/25.

The 25 templates I_k depend only on geometry (the cell-center positions ξ_k), not on the densities, so they are computed once per run and reused for every element. With uniform density this 5 by 5 midpoint rule reproduces the analytical Q4 stiffness to the order of the midpoint rule, which is high accuracy here because the integrand is a low-order polynomial. The key point: the dominant approximation in MTOP is the coarse displacement mesh, not this quadrature. The quadrature error is second order by comparison.

### Q3. Why does the native MTOP compliance (219) differ from the fine-mesh re-analysis (224.7)? Which number should be compared?

They are the same density field evaluated by two different analysis models. The native MTOP value (219.069) uses the coarse 120 by 40 displacement model with the midpoint-integrated stiffness; the fine value (224.696) uses the standard 600 by 200 FE model. A displacement-based FE model is always too stiff, because the assumed displacement field is only a subspace of all admissible deformations, so the structure cannot deform as freely as it really would. A coarser mesh is a smaller subspace, hence stiffer, hence lower compliance. The roughly 2.5 % gap is therefore a discretization effect, not an optimization error.

The fine-mesh re-analysis is the correct quantity for comparing mechanical quality, because it evaluates every design (classical, MTOP, coarse baseline) with one common, converged model. The native compliance is only an internal diagnostic of the coarse analysis. That is why every headline comparison in the report uses fine-mesh compliance.

### Q4. MTOP can produce non-physical density patterns inside an element. Why does that happen, and what prevents it here?

This is the known weakness of multiresolution schemes. The coarse analysis only sees the 25 densities through their integrated contribution to the 8 by 8 element stiffness. Many different fine-scale density patterns give nearly the same K_e, and some can look artificially efficient: a pattern aligned with the coarse element's strain modes can be credited with high stiffness even though, re-analyzed on a fine mesh, it is a disconnected or checkerboard-like layout. The coarse displacement field simply cannot resolve the difference.

Two things prevent it here. First, the cone filter (radius 24 density cells, about 5 coarse elements) correlates neighboring densities and enforces a minimum length scale much larger than one coarse element, so high-frequency intra-element patterns are filtered out before they can be exploited. Second, the fine-mesh re-analysis is the safety net: if MTOP had exploited an artificial pattern, the fine-mesh compliance would jump. It does not (0.045 %), which confirms the filter is doing its job.

### Q5. How is the per-cell sensitivity of compliance derived, and how does it differ from the classical element sensitivity?

Compliance is self-adjoint (see Q7), so the general result is dC/dρ = − u^T (dK/dρ) u. Only the cell's own template contributes to dK/dρ_ek, because K_e = Σ_k [...] I_k and the cells are independent: dK_e/dρ_ek = p (E0 − E_min) ρ_ek^(p−1) I_k. Therefore

dC/dρ_ek = − p (E0 − E_min) ρ_ek^(p−1) u_e^T I_k u_e.

The classical element sensitivity is dC/dρ_e = − p (E0 − E_min) ρ_e^(p−1) u_e^T K0 u_e, a single strain-energy density per element. The MTOP version replaces the element strain energy u_e^T K0 u_e by a per-cell strain energy u_e^T I_k u_e evaluated at the cell center. So each of the 25 cells gets its own gradient even though they share the same element displacement vector u_e. It is a clean chain-rule result, which we verified by finite differences (Taylor slope 2.00).

### Q6. Why is the measured speedup 21 times and not the 40 times of the FE solve itself?

The phase breakdown answers this. In the classical OC run the FE solve (factorizing and solving K u = f on about 240,000 DOF) is 82 % of the time, and MTOP shrinks exactly that: the coarse solve on about 9,900 DOF is roughly 40 times faster (about 215 s down to about 5 s). But the other phases do not shrink. The cone filter, the assembly of K_e from 25 templates, and the OC update all still operate on the full 600 by 200 density mesh, and their cost is essentially identical in both formulations.

This is Amdahl's law: once you make the dominant 82 % part 40 times cheaper, the remaining 18 % becomes the new bottleneck and caps the total speedup at about 21. With density filtering the fixed per-density-cell work is even larger (an extra chain-rule convolution and a filtered-volume check), so the speedup falls further, to about 5.

## B. Sensitivity analysis

### Q7. Derive the compliance sensitivity. Why is the problem self-adjoint, and what does that save?

Compliance is C = f^T u with state equation K(ρ) u = f. Differentiating the state equation, and using that the load is design-independent (df/dρ = 0):

(dK/dρ) u + K (du/dρ) = 0, so du/dρ = − K^(−1) (dK/dρ) u.

Then dC/dρ = f^T (du/dρ) = − f^T K^(−1) (dK/dρ) u = − (K^(−1) f)^T (dK/dρ) u = − u^T (dK/dρ) u, using K symmetric and K^(−1) f = u.

The general adjoint method introduces an adjoint vector λ solving K^T λ = dC/du. Here the objective is C = f^T u, so dC/du = f, and K is symmetric, so the adjoint equation is K λ = f, which is identical to the state equation: λ = u. That is what self-adjoint means: the adjoint solution is the displacement itself, so no extra linear solve is needed for the gradient. With SIMP, dK/dρ_e = p ρ_e^(p−1) (E0 − E_min) K0, giving the familiar form in Q5. If the objective were instead a displacement at a single point, the problem would not be self-adjoint and a separate adjoint solve would be required.

### Q8. Why does the Taylor remainder test have slope 2, and why is it stronger than comparing the gradient to a forward difference?

For a perturbation h d in a direction d, a Taylor expansion gives C(x + h d) = C(x) + h (gradC . d) + O(h^2). The zeroth-order remainder, |C(x + h d) − C(x)|, is O(h): it shrinks linearly. The first-order remainder, |C(x + h d) − C(x) − h (gradC . d)|, subtracts off the term that uses our analytical gradient. If the gradient is correct, what remains is O(h^2), which is a slope of 2 on a log-log plot. If the gradient were wrong, the h (gradC . d) term would fail to cancel the true linear term and the remainder would stay O(h), slope 1.

That is why it is stronger than a forward-difference comparison. A forward difference is only a first-order accurate estimate of the gradient, so agreement to a few digits is the best you can hope for, and a sign or scaling error can hide. The Taylor test checks the convergence rate: obtaining slope 2 (we obtain 2.00) confirms not just the magnitude but that the analytical gradient is exactly the first-order term. It is also direction-based, so a single sweep validates the whole gradient vector at once.

### Q9. Why do the finite-difference errors get worse again at very small step sizes?

A finite difference has two competing error sources. Truncation error falls as h decreases (O(h) for a forward difference, O(h^2) for a central difference). Round-off error grows as h decreases: the numerator C(x + h) − C(x) becomes the subtraction of two nearly equal numbers, and in double precision you lose significant digits through cancellation; that error scales like machine epsilon divided by h. The total error is therefore U-shaped on a log-log plot: truncation-dominated for large h, round-off-dominated for small h, with an optimal step size in between (around 1e-5 to 1e-6 for our checks, where the relative error bottoms out near 1e-6). Seeing this U-shape is itself a sanity check: if the curve did not turn up at small h, it would suggest the analytical gradient is merely tracking the finite difference rather than being an independent quantity.

## C. SIMP and filtering

### Q10. Why penalization power p = 3? What happens at p = 1 or at large p, and why does SIMP need a filter?

SIMP sets E(ρ) = E_min + ρ^p (E0 − E_min). The purpose of p > 1 is to make intermediate density uneconomical. At ρ = 0.5 with p = 3, the normalized stiffness is only E/E0 = 0.125: half the material buys one-eighth of the stiffness. Volume scales linearly with ρ but stiffness scales as ρ^p, so for the optimizer gray material is a bad deal, and the design is driven toward 0 and 1.

At p = 1 there is no penalization: the stiffness-to-volume ratio is constant, gray is as efficient as solid, and the optimum is a blurry gray field rather than a discrete structure. At very large p the 0/1 contrast is sharper but the problem becomes strongly non-convex with many poor local minima, and convergence suffers. p = 3 is the established compromise.

SIMP still needs a filter because the discretized problem is otherwise ill-posed. Without a filter the optimizer exploits the FE discretization to form checkerboard patterns, which are artificially stiff, and the solution is mesh-dependent: a finer mesh just produces finer members and the design never converges under mesh refinement. The filter imposes a minimum length scale, which restores mesh-independence and removes checkerboards.

### Q11. Sensitivity filtering versus density filtering: what is the difference, and why does sensitivity filtering have no clean physical density?

Both use the same cone weights H_ij = max(0, r_min − distance between cells i and j).

Density filtering filters the design variables: the physical density is ρ_i = (Σ_j H_ij x_j) / (Σ_j H_ij). This is a genuine, well-defined field. The analysis and the volume constraint use ρ, and the chain rule dC/dx = (dρ/dx)^T dC/dρ is exact. The filter is part of the objective.

Sensitivity filtering keeps the densities equal to the design variables (ρ_i = x_i) and instead smooths the gradient: the filtered sensitivity is a weighted average of the neighboring x_j (dC/dx_j), divided by max(γ, x_i). It is a heuristic that modifies the search direction, not the objective. There is no functional whose exact gradient is the filtered sensitivity, so there is no physical density distinct from x. It is cheap and works well, but because it is not a consistent gradient it sits less comfortably inside a formal optimizer such as MMA, which is part of why the MMA sensitivity-filtered run needs more iterations than its OC counterpart.

### Q12. Why does the cone filter impose a minimum length scale, and what does it protect against?

The cone kernel replaces each cell value by a weighted average over a disc of radius r_min (24 density cells here), with weights decaying linearly to zero at the radius. Any density variation with a wavelength shorter than about r_min is averaged out: the filtered field cannot contain features finer than the filter radius. That is the minimum length scale.

It protects against the two pathologies of unfiltered SIMP. Checkerboarding: alternating solid and void on adjacent elements is artificially stiff because of the bilinear FE displacement field, and it has a one-element wavelength, so the filter destroys it. Mesh dependence: without a filter, refining the mesh lets the optimizer insert ever-finer members, so the design never settles under mesh refinement; fixing a physical r_min decouples the feature size from the mesh. As a bonus it gives a manufacturability guarantee, a minimum member width.

### Q13. The OC plus density-filter runs never meet the design-change tolerance. Explain precisely why, and is the design actually converged?

The stopping test is max_i |Δx_i| < 0.01. The OC update is multiplicative and per-variable: each x_i is independently rescaled by the square root of B_i/λ. The density filter, however, couples the variables: every design-variable change is spread over a neighborhood before it reaches the physical density. So the optimizer nudges a variable, the filter smears that change, and at member boundaries the per-variable OC rescaling and the filter smoothing work against each other. The result is a small, persistent high-frequency oscillation of the design variables at the boundaries: max |Δx| stalls around 0.04 to 0.07 and never drops below 0.01.

This is not a non-converged design. The compliance is stable to well within 0.1 % of its final value, which is why we add a compliance-plateau criterion: terminate when the relative spread of the compliance over a 20-iteration window falls below 1e-4. The structure is converged; only the design-variable vector keeps jittering. It is an artifact of pairing a heuristic multiplicative update with a consistent filter, and switching to MMA does not fully cure it either.

## D. Optimizers

### Q14. Explain the OC update rule. Where does the multiplicative form come from, and what does the bisection on the Lagrange multiplier do?

For minimum compliance with a single volume constraint, the optimality (KKT) condition is that, at every interior variable, the ratio B_i = − (filtered dC/dx_i) / (dV/dx_i) equals the Lagrange multiplier λ of the volume constraint. Physically B_i is the stiffness gained per unit volume at cell i; at the optimum this is equalized across all interior cells.

The OC update is a fixed-point iteration toward B_i = λ: x_i^new = x_i (B_i/λ)^(1/2), clipped to the move-limited interval [x_i − m, x_i + m] and to [0, 1], with move limit m = 0.2. If a cell is under-used (B_i > λ) it grows; if over-used it shrinks; at the optimum it is stationary. The exponent one half (the square root) damps the step so the iteration is stable.

λ is not known in advance; it must make the volume constraint active. For any trial λ the update produces a design with some volume, and volume decreases monotonically with λ, so we bisect on λ until the volume equals V* = 0.5. The constraint is always active for compliance minimization, because more material never reduces stiffness. So each OC iteration is an inner bisection wrapped around the multiplicative update.

### Q15. What does MMA approximate, and why are the asymptotes adaptive? Why is MMA better suited to Heaviside projection than OC?

MMA (the method of moving asymptotes, due to Svanberg) is a sequential convex programming method. At each iterate x_k it replaces the objective and every constraint by a separable, convex approximation that matches the function value and gradient at x_k. Each term has the form a/(U_k − x) + b/(x − L_k), with moving asymptotes L_k < x_k < U_k. This is convex and diverges to plus infinity at the asymptotes, so the asymptotes act as a trust region. The convex subproblem is solved exactly by an interior-point method, giving x_{k+1}.

The asymptotes are adaptive because they encode the local curvature and trust. If a variable oscillates between iterations the asymptotes are pulled inward, making the next step more conservative; if it moves monotonically they are pushed outward, allowing a larger step.

MMA suits Heaviside projection because the projection makes the objective strongly non-linear in x (a steep tanh), and MMA is a general gradient-based method that handles arbitrary differentiable objectives and constraints. OC, in contrast, is a heuristic derived specifically for the plain compliance and volume structure; its multiplicative form has no mechanism for the sharp projection non-linearity. OC also implicitly assumes a consistent gradient, whereas the sensitivity filter supplies only a heuristic one.

### Q16. Why does the Heaviside-projected design have lower compliance than the gray density-filtered design at the same volume fraction?

Because gray (intermediate-density) material is penalized by SIMP and is therefore mechanically inefficient. A density-filtered design has wide gray transition bands at every member edge, since by construction the filter cannot make a boundary sharper than r_min. At ρ = 0.5 the SIMP stiffness is only 0.125 of E0, so that boundary material contributes volume but very little stiffness.

Heaviside projection pushes those transition cells toward 0 or 1. The volume fraction is held at 0.5 by the constraint, so the same amount of material is redistributed out of inefficient gray bands and into either solid members (where it carries load at the full E0) or void. Concentrating the material where it actually works lowers the compliance, by about 17 % here (240 down to 199). The projected design is also closer to a true black-and-white structure, which is what the discrete topology-optimization problem really wants; the gray design is partly an artifact of the continuous relaxation.

### Q17. Why the beta-continuation? What would go wrong if you started at beta = 16?

Beta controls the sharpness of the Heaviside tanh projection: small beta is almost linear (barely any projection), large beta approaches a true step. The projection derivative dρ/dρ_tilde grows large and very localized as beta increases: it is near zero everywhere except in a thin band around the threshold η.

If you start at beta = 16, the optimizer faces an almost discontinuous objective from the first iteration. Gradients vanish over most of the domain, the convex MMA subproblem is a poor local model, and the run locks into a bad local minimum close to the still-gray starting design. Continuation, beta starting at 1 and doubling every 50 iterations up to 16, lets the optimizer first find a good gray layout while the problem is smooth and nearly convex, then gradually sharpen the boundaries while staying in the basin of that good design. It is a homotopy or numerical-continuation strategy: solve an easy problem and deform it slowly into the hard one.

## E. Results, interpretation, and limitations

### Q18. Your plain coarse-resolution baseline reaches within 0.6 % of the fine design. Doesn't that undermine the case for MTOP?

It is a fair question, and we added that baseline deliberately to be honest about it. For this problem, 2D with a fairly large filter radius (24 density cells), the minimum feature size is large, so even a plain 120 by 40 model resolves the load path well: re-analyzed on the fine mesh it is only 0.6 % behind.

But the comparison is not like for like. That coarse design is locked to 120 by 40 resolution and is measurably grayer (grayness index 0.265 against a sharper MTOP design): part of its low compliance is unpenalized intermediate material, and it cannot be post-processed or manufactured at fine resolution. MTOP keeps the full 600 by 200 design description, a genuinely fine and near-black-and-white layout, while paying only the coarse solve. So MTOP is not competing with the coarse model on compliance; it delivers a fine-resolution design at coarse-solve cost. The case for MTOP gets stronger exactly where the coarse baseline gets weaker: smaller features, a smaller filter radius, and 3D (see Q19).

### Q19. Is the 21 times speedup hardware-independent? How would it change in 3D?

No, the absolute timings and the exact factor are run-specific: they depend on the machine, the MATLAB release, the sparse solver, and the linear-algebra libraries. The report is explicit that the speedup is a run-specific comparison made under identical conditions, not a hardware constant; both formulations were timed back to back on the same machine and release.

The trend, though, is robust and would grow in 3D. The speedup comes from shrinking the FE solve. For a sparse direct solver, the factorization cost scales roughly as O(n^1.5) in 2D but about O(n^2) in 3D, because the fill-in is far worse. So the solve dominates the iteration even more heavily in 3D than the 82 % we measure here, and reducing the analysis mesh by the same factor buys a larger share of the total. The operations that cap our speedup (filter, assembly, optimizer on the fine density mesh) are comparatively cheap next to a 3D solve. So 21 times is a conservative, 2D figure.

### Q20. Beyond the finite-difference checks, how do you know the MTOP implementation is correct?

There are four independent layers of evidence. First, sensitivity verification: three finite-difference strategies on three derivative chains, with a Taylor-remainder slope of 2.00, which validates the gradients and therefore that the optimizer descends the true objective. Second, consistency limits: with uniform density the 25-template MTOP stiffness reproduces the analytical Q4 stiffness, and the volume sensitivity is recovered to machine precision. Third, cross-model agreement: the MTOP design re-analyzed on the independent classical 600 by 200 FE model gives a compliance within 0.045 % of the classical design, so two separately coded analysis paths agree. Fourth, physical plausibility: the optimized layout is the textbook MBB truss (a top chord in compression, a bottom chord in tension, a triangulated web), the volume constraint is satisfied exactly, and the convergence histories behave as expected. Any indexing or chain-rule bug would break at least one of these checks; in fact the component-wise finite-difference check, which localizes an error to a single design variable, is what we used to find bugs during development.

## Key numbers to have ready

Optimized-design compliance (fine-mesh re-analysis unless noted) and cost:

- Classical OC, sensitivity filter: C = 224.6, 94 iterations, 263 s
- MTOP OC, sensitivity filter: C = 224.7 (+0.045 %), 95 iterations, 12 s, speedup 21.4 times
- Classical OC, density filter: C = 241.1, 237 iterations, 693 s
- MTOP OC, density filter: C = 241.0, 299 iterations, 137 s, speedup 5.07 times
- MTOP MMA, sensitivity filter: C = 226.5, 250 iterations, 187 s
- MTOP MMA, density filter: C = 240.1, 419 iterations, 254 s
- MTOP MMA, Heaviside: C = 200.0 / 199.3 / 200.6 for the thresholds 0.3 / 0.5 / 0.7
- Coarse-resolution baseline (120 by 40 analysis equal to density mesh): C = 226.0, about 0.6 % above the MTOP design

Method and verification:

- Free degrees of freedom: classical about 241,000, MTOP analysis about 9,900 (mesh-size factor about 24)
- FE solve: about 215 s classical against about 5 s MTOP, roughly 40 times faster
- FE solve as a share of classical OC time: 82 %
- SIMP penalization p = 3; at the volume fraction ρ = 0.5 the normalized stiffness is E/E0 = 0.125
- Cone filter radius r_min = 24 density cells; target volume fraction V* = 0.5; OC move limit m = 0.2
- Taylor-remainder slope about 2.00; best finite-difference relative error about 1e-6
- Heaviside continuation: beta from 1 to 16, doubled every 50 iterations
