function dof = dof_solid8(Nodenumbers)
%DOF_SOLID8   Element degrees of freedom for a solid8 element.
%
%   dof = dof_solid8(NodeNum) builds the vector with the 
%   labels of the degrees of freedom for which stiffness is present in the 
%   solid8 element. 
%
%   NodeNum Node definitions           [NodID1 NodID2 ... NodIDn]   (1 * 4)
%   dof         Degrees of freedom                                 (24 * 1)       
%
%   See also GETDOF.

if nargin<1
   Nodenumbers = (1:8); 
end

dof=zeros(24,1);
dof(1:3)  = Nodenumbers(1) + [0.01 0.02 0.03]';
dof(4:6)  = Nodenumbers(2) + [0.01 0.02 0.03]';
dof(7:9)  = Nodenumbers(3) + [0.01 0.02 0.03]';
dof(10:12)= Nodenumbers(4) + [0.01 0.02 0.03]';
dof(13:15)= Nodenumbers(5) + [0.01 0.02 0.03]';
dof(16:18)= Nodenumbers(6) + [0.01 0.02 0.03]';
dof(19:21)= Nodenumbers(7) + [0.01 0.02 0.03]';
dof(22:24)= Nodenumbers(8) + [0.01 0.02 0.03]';
