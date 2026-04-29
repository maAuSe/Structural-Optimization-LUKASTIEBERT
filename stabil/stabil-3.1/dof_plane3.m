function dof = dof_plane3(NodeNum)

%DOF_PLANE3   Element degrees of freedom for a plane3 element.
%
%   dof = dof_plane3(NodeNum) builds the vector with the 
%   labels of the degrees of freedom for which stiffness is present in the 
%   plane3 element.
%
%   NodeNum Node definitions           [NodID1 NodID2 NodID3]   (1 * 3)
%   dof     Degrees of freedom                                  (6 * 1)       
%
%   See also GETDOF.

nNode=3;

if nargin<1
    NodeNum=1:nNode;
end

dof=zeros(6,1);
dof(1:2)=NodeNum(1)+[0.01 0.02].';
dof(3:4)=NodeNum(2)+[0.01 0.02].';
dof(5:6)=NodeNum(3)+[0.01 0.02].';
end

