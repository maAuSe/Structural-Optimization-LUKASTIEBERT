# Examination assignment

## Multiresolution topology optimization (MTOP)

Mattias Schevenels

## General guidelines

The evaluation of the course *Structural Optimization* is based on a group assignment (groups of three students). Students are required to submit a written report and give an oral presentation.

### Dates and deadlines

- Start of the assignment: 3 April 2026.
- Report submission deadline: 15 May 2026.
- Presentations: 22 May 2026.

### Evaluation

- The evaluation will consider the quality of the methodology, correctness of the results, clarity of the report and presentation, and the demonstrated understanding of the work.
- In addition, students are expected to demonstrate an understanding of the course material relevant to their assignment. During the oral defense, questions may address theoretical concepts and methods that are relevant to or could be meaningfully applied in the context of their assignment.

## Report guidelines

- The report can be written in Dutch or in English.
- The report should be concise but complete. All relevant information must be included, in particular:
  - Problem statement.
  - Definition of the structural model: geometry, boundary conditions, loading.
  - Formulation of the optimization problem: design variables, objective function, constraints.
  - Description of the solution technique including relevant parameter values. Use mathematical symbols and equations; do not include MATLAB code.
  - Verification of the implementation, in particular the verification of the sensitivities (e.g. using finite differences).
  - Figures of the results: optimized designs, convergence histories.
  - Physical interpretation of the results: explanation of intermediate and final results from a physical perspective, demonstration of the influence of certain parameters, comparison of results obtained using different or similar optimization methods, etc.
  - Transparency statement on the use of GenAI (see below).
  - Bibliography containing all consulted references.
- The report should be submitted by email (mattias.schevenels@kuleuven.be) in PDF format, together with a ZIP archive containing all MATLAB files and a readme.txt file with a brief overview of the MATLAB files. The MATLAB files are requested to ensure reproducibility of the results: they should generate all figures included in the report. As part of the assessment, the MATLAB files may be executed to verify if this is the case. The MATLAB code itself will, in principle, not be assessed. Therefore, all information required to understand the methodology and results must be included in the report.

## Oral presentation and defense guidelines

- Every group will give a presentation of 10 minutes, followed by 10 minutes for questions.
- The presentation is given in English.
- The slides should be submitted by email before the start of the presentations.
- The target audience consists of the other students who followed the structural optimization course: it is not necessary to re-explain the basics of structural optimization, but it is important to explain the problem statement and specific methods or techniques used in the context of your assignment.
- Discuss the solution strategy. Do not show MATLAB code, limit the use of equations, and try to convey the message as much as possible by means of figures showing intermediate or final results.
- Focus on the physical interpretation of the results.
- Briefly explain how you used GenAI and how you ensured that you fully understood and validated the results.
- All group members are expected to contribute and should be able to explain all parts of the work.
- All students are expected to attend all presentations.

## GenAI guidelines

GenAI has become a powerful tool, not only for programming tasks but also for numerical modeling, and therefore also in the context of this course. The use of GenAI for this assignment is allowed. However, it must be used thoughtfully and appropriately:

- From an engineering perspective: you remain responsible for the results. You must be able to explain what is being computed, how it is computed, and why this is the correct approach, and you have to ensure that the results are correct.
- From an educational perspective: you must demonstrate that you have achieved the learning objectives. In the report, during the presentation, and in the discussion afterwards, you must show that you are in control of the work.
- From a scientific integrity perspective: you must be transparent about how you have used GenAI.

These three aspects are closely related and all require a critical attitude which you must also demonstrate.

The use of GenAI is also allowed for writing the report and preparing the presentation. Transparency remains essential here as well.

## Guidance and feedback

- For guidance, you can make an appointment by email (mattias.schevenels@kuleuven.be).
- Every group is expected to request guidance at least once, to discuss the results before including them in the report.
- After the announcement of the results, it is possible to make an appointment (by email) for feedback.

## Problem statement

This assignment focuses on the so-called multiresolution topology (MTOP) optimization approach introduced by Nguyen et al. [1]. The basic idea of the MTOP approach is to represent the design by means of a density mesh with a fine resolution, and to perform the analysis on a coarser finite element mesh. This leads to a smaller system of equilibrium equations, which should result in lower computation times. However, the original paper [1] does not report computation times or speedup factors. The aim of this assignment is to assess the computational efficiency of the MTOP approach for 2D minimum compliance topology optimization.

## Assignment

1. Start from the MBB beam code from the lectures. Solve the problem with the optimality criteria method, for a mesh of 600 × 200 elements, a volume fraction of 0.5, a penalization power of 3, and sensitivity filtering with a filter radius of 24 elements. Plot the convergence history (evolution of the compliance and the volume fraction) in terms of (1) the iteration number and (2) the computation time.
2. Solve the same problem with the MTOP approach. Use a finite element mesh of 120 × 40 elements and a density mesh of 5 × 5 elements per finite element. Visualize the difference between the resulting design and the design obtained in step 1. Plot the convergence history and compare it with the results from step 1.
3. Investigate how the MTOP approach performs as compared to the classical approach when using density filtering.
4. Switch to MMA and investigate the performance of the MTOP approach when using sensitivity filtering, density filtering, and Heaviside projection (using different values for the projection threshold).

## References

[1] T. Nguyen, G.H. Paulino, J. Song, and H.L. Chua. A computational paradigm for multiresolution topology optimization (MTOP). *Structural and Multidisciplinary Optimization*, 41(4):525–539, 2010.