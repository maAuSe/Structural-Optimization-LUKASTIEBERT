function dof = dof_shell6(NodeNum)

%DOF_SHELL6   Element degrees of freedom for a shell6 element.
%
%   dof = dof_shell6(NodeNum) builds the vector with the 
%   labels of the degrees of freedom for which stiffness is present in the 
%   shell6 element.
%
%   NodeNum Node definitions           [NodID1 NodID2 ... NodIDn]   (1 * 6)
%   dof     Degrees of freedom                                     (36 * 1)       
%
%   See also GETDOF.

% Miche Jansen
% 2009

if nargin<1
NodeNum = (1:6);
end

NodeNum = NodeNum(:).';
dof2=(0.01:0.01:0.06).';
dof = NodeNum(ones(6,1),1:6)+dof2(:,ones(1,6));
dof = dof(:);
end



