function [UeGCS,L] = elemdisp(Type,NodeNum,DOF,U)

%ELEMDISP  Select the element displacements from the global displacement vector.
%
%   UeGCS = elemdisp(Type,NodeNum,DOF,U) selects the element
%   displacements from the global displacement vector.
%
%   Type       Element type e.g. 'beam','truss', ...
%   NodeNum    Node numbers (1 * nNodes)
%   DOF        Degrees of freedom  (nDOF * 1)
%   U          Displacements (nDOF * nLC)
%   UeGCS      Element displacements
%
%   See also ELEMFORCES, DOF_BEAM, DOF_TRUSS.

% David Dooms
% September 2008

dofelem=eval(['dof_' Type '(NodeNum)']);

L=selectdof(DOF,dofelem);

UeGCS=L*U;       
