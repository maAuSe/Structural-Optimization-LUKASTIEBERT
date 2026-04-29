function dof = dof_shell2(NodeNum)

%DOF_SHELL2   Element degrees of freedom for a SHELL2 element.
%
%   DOF = dof_shell2(NodeNum) returns a vector with the degrees of freedom for
%   a SHELL2 element.
%
%   NodeNum Node definitions           [NodID1 NodID2]  (1 * 2)
%   DOF     Degrees of freedom                          (12 * 1)

% Mattias Schevenels
% April 2020

if nargin<1
    NodeNum = (1:2);
end

dof=zeros(6,1);
dof(1:3,1)=NodeNum(1)+[0.01 0.02 0.06].';
dof(4:6,1)=NodeNum(2)+[0.01 0.02 0.06].';
end
