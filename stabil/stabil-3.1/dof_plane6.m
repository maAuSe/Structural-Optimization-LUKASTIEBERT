function dof = dof_plane6(NodeNum)

%DOF_PLANE6   Element degrees of freedom for a plane6 element.
%
%   dof = dof_plane6(NodeNum) builds the vector with the 
%   labels of the degrees of freedom for which stiffness is present in the 
%   plane6 element.
%
%   NodeNum Node definitions           [NodID1 NodID2 ... NodIDn]   (1 * 6)
%   dof     Degrees of freedom                                      (12 * 1)       
%
%   See also GETDOF.

% Jef Wambacq
% 2017

nNode=6;

if nargin<1
    NodeNum=1:nNode;
end

NodeNum = NodeNum(:).';
dof2 = [0.01;0.02];
dof = [NodeNum(1:nNode);NodeNum(1:nNode)] + dof2(:,ones(nNode,1));
dof = dof(:);

end



