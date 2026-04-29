function dof = dof_beam(NodeNum)

%DOF_BEAM   Element degrees of freedom for a beam element.
%
%   DOF = dof_beam(NodeNum) builds the vector with the 
%   builds the vector with the degrees of freedom for the beam element.
%
%   NodeNum Node definitions           [NodID1 NodID2]  (1 * 2)
%   DOF     Degrees of freedom                          (12 * 1)
%
%   See also GETDOF.

% David Dooms
% March 2008

if nargin<1
    NodeNum = (1:2);
end

dof=zeros(12,1);
dof(1:6,1)=NodeNum(1)+[0.01:0.01:0.06].';
dof(7:12,1)=NodeNum(2)+[0.01:0.01:0.06].';
end
