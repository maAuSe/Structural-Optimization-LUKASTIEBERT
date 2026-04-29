function dof = dof_shell8(NodeNum)

%DOF_SHELL8   Element degrees of freedom for a shell8 element.
%
%   dof = dof_shell8(NodeNum) builds the vector with the 
%   labels of the degrees of freedom for which stiffness is present in the 
%   shell8 element.
%
%   NodeNum Node definitions           [NodID1 NodID2 ... NodIDn]   (1 * 8)
%   dof     Degrees of freedom                                     (48 * 1)       
%
%   See also GETDOF.

% Miche Jansen
% 2009

if nargin<1
NodeNum = (1:8);
end

NodeNum = NodeNum(:).';
dof2=(0.01:0.01:0.06).';
dof = NodeNum(ones(6,1),1:8)+dof2(:,ones(1,8));
dof = dof(:);
end



