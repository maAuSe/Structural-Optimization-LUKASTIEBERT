function dof = dof_mass(NodeNum)
%DOF_MASS   Element degrees of freedom for a mass element.
%
%   dof = dof_truss(NodeNum) builds the vector with the 
%   labels of the degrees of freedom for which stiffness is present in the 
%   mass element.
%
%   NodeNum Node numbers           [NodID1] (1)
%   dof         Degrees of freedom  (6 * 1)
%
%   See also GETDOF.

dof=zeros(6,1);
dof(1:6,1)=NodeNum(1)+[0.01; 0.02; 0.03; 0.04; 0.05; 0.06];
