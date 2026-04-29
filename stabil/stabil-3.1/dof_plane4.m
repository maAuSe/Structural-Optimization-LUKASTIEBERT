function dof = dof_plane4(NodeNum)

%DOF_PLANE4   Element degrees of freedom for a plane4 element.
%
%   dof = dof_plane4(NodeNum) builds the vector with the 
%   labels of the degrees of freedom for which stiffness is present in the 
%   plane4 element.
%
%   NodeNum Node definitions           [NodID1 NodID2 ... NodIDn]   (1 * 4)
%   dof     Degrees of freedom                                      (8 * 1)       
%
%   See also GETDOF.

% Miche Jansen
% 2012

if nargin<1
NodeNum = (1:4);
end

NodeNum = NodeNum(:).';
dof2 = [0.01;0.02];
dof = [NodeNum(1:4);NodeNum(1:4)] + dof2(:,[1,1,1,1]);
dof = dof(:);

end