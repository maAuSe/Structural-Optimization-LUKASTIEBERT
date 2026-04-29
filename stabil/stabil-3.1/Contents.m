% TOOLBOX
%
% Files
%   getdof          - Get the vector with the degrees of freedom of the model.
%   asmkm           - Assemble stiffness and mass matrix.
%   removedof       - Remove DOF with Dirichlet boundary conditions equal to zero.
%   addconstr       - Add constraint equations to the stiffness matrix and load vector.
%   tconstr         - Return matrices to apply constraint equations.
%   nodalvalues     - Construct a vector with the values at the selected DOF.
%   elemloads       - Equivalent nodal forces.
%   accel           - Compute the distributed loads due to an acceleration.
%   elemforces      - Compute the element forces.
%   elemdisp        - Select the element displacements from the global displacement vector.
%   selectdof       - Select degrees of freedom.
%   unselectdof     - Unselect degrees of freedom.
%   selectnode      - Select nodes by location.
%   reprow          - Replicate rows from a matrix.
%   multdloads      - Combine distributed loads.

%   dof_beam        - Element degrees of freedom for a beam element.
%   ke_beam         - Beam element stiffness and mass matrix in global coordinate system.
%   kelcs_beam      - Beam element stiffness and mass matrix in local coordinate system.
%   trans_beam      - Transform coordinate system for a beam element.
%   loads_beam      - Equivalent nodal forces for a beam element in the GCS.
%   loadslcs_beam   - Equivalent nodal forces for a beam element in the LCS.
%   accel_beam      - Compute the distributed loads for a beam due to an acceleration.
%   forces_beam     - Compute the element forces for a beam element.
%   forceslcs_beam  - Compute the element forces for a beam element in the LCS.
%   nelcs_beam      - Shape functions for a beam element.
%   nedloadlcs_beam - Shape functions for a distributed load on a beam element.
%   coord_beam      - Coordinates of the beam elements for plotting.
%   disp_beam       - Return matrices to compute the displacements of the deformed beams.
%   dispgcs2lcs_beam- Transform the element displacements to the LCS for a beam.
%   fdiagrgcs_beam  - Return matrices to plot the forces in a beam element.
%   fdiagrlcs_beam  - Force diagram for a beam element in LCS.
%   sdiagrgcs_beam  - Return matrices to plot the stresses in a beam element.
%   sdiagrlcs_beam  - Stress diagram for a beam element in LCS.


%   dof_truss       - Element degrees of freedom for a truss element.
%   ke_truss        - Truss element stiffness and mass matrix in global coordinate system.
%   kelcs_truss     - Truss element stiffness and mass matrix in local coordinate system.
%   trans_truss     - Transform coordinate system for a truss element.
%   loads_truss     - Equivalent nodal forces for a truss element in the GCS.
%   accel_truss     - Compute the distributed loads for a truss due to an acceleration.
%   forces_truss    - Compute the element forces for a truss element.
%   forceslcs_truss - Compute the element forces for a truss element in the LCS.
%   coord_truss     - Coordinates of the truss elements for plotting.
%   disp_truss      - Return matrices to compute the displacements of the deformed trusses.
%   dispgcs2lcs_truss- Transform the element displacements to the LCS for a truss.
%   fdiagrgcs_truss - Return matrices to plot the forces in a truss element.
%   sdiagrgcs_truss - Return matrices to plot the stresses in a truss element.

%   plotnodes       - Plot the nodes.
%   plotelem        - Plot the elements.
%   plotdisp        - Plot the displacements.
%   plotforc        - Plot the forces.
%   plotstress      - Plot the stresses.
%   animdisp        - Animate the displacements.
%   getmovie        - Get the movie from a figure where an animation has been played.     
%   printdisp       - Display the displacements in the command window.
%   printforc       - Display the forces in the command window.

%   eigfem          - Compute the eigenmodes and eigenfrequencies of the finite element model.
%   msupf           - Modal superposition in the frequency domain.
%   msupt           - Modal superposition in the time domain.
%   cdiff           - Direct time integration for dynamic systems - central diff. method.
%   newmark         - Direct time integration for dynamic systems - Newmark method
%   wilson          - Direct time integration for dynamic systems - Wilson-theta method
