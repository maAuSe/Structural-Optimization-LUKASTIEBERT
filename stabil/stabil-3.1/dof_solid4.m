function dof = dof_solid4(Nodenumbers)
%DOF_SOLID4   Element degrees of freedom for a solid4 element.
%
%   dof = dof_solid4(NodeNum) builds the vector with the 
%   labels of the degrees of freedom for which stiffness is present in the 
%   solid4 element. 
%
%   NodeNum   Node definitions    [NodID1 NodID2 NodID3 NodID4]   (1 * 4).
%   dof       Degrees of freedom                                 (12 * 1).      
%
%   See also GETDOF.

dof=zeros(12,1);
dof(1:3)  = Nodenumbers(1)  + [0.01 0.02 0.03]';
dof(4:6)  = Nodenumbers(2)  + [0.01 0.02 0.03]';
dof(7:9)  = Nodenumbers(3)  + [0.01 0.02 0.03]';
dof(10:12)= Nodenumbers(4)  + [0.01 0.02 0.03]';
