function dof = dof_truss(NodeNum)

%DOF_TRUSS   Element degrees of freedom for a truss element.
%
%   dof = dof_truss(NodeNum) builds the vector with the degrees of freedom 
%   for the truss element.
%
%   NodeNum Node numbers           [NodID1 NodID2] (1 * 2)
%   dof         Degrees of freedom  (6 * 1)
%
%   See also GETDOF.

% David Dooms
% March 2008

if nargin<1
NodeNum = (1:2);
end

dof=zeros(6,1);
dof(1:3,1)=NodeNum(1)+[0.01; 0.02; 0.03];
dof(4:6,1)=NodeNum(2)+[0.01; 0.02; 0.03];
