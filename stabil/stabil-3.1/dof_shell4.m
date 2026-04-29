function dof = dof_shell4(NodeNum)

%DOF_SHELL4   Element degrees of freedom for a shell4 element.
%
%   dof = dof_shell4(NodeNum) builds the vector with the 
%   labels of the degrees of freedom for which stiffness is present in the 
%   shell4 element.
%
%   NodeNum Node definitions           [NodID1 NodID2 ... NodIDn]   (1 * 4)
%   dof         Degrees of freedom                                 (24 * 1)       
%
%   See also GETDOF.

% Miche Jansen
% 2009

if nargin<1
NodeNum = (1:4);
end

NodeNum = NodeNum(:).';
dof2 = (0.01:0.01:0.06).';
dof = NodeNum([1 1 1 1 1 1],1:4) + dof2(:,[1 1 1 1]);
dof = dof(:);

end