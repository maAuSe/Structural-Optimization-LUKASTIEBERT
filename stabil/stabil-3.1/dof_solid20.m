function dof = dof_solid20(Nodenumbers)
%DOF_SOLID20   Element degrees of freedom for a solid20 element.
%
%   dof = dof_solid20(NodeNum) builds the vector with the 
%   labels of the degrees of freedom for which stiffness is present in the 
%   solid20 element. 
%
%   NodeNum Node definitions           [NodID1 NodID2 ... NodIDn]   (1 * 4)
%   dof         Degrees of freedom                                 (60 * 1)       
%
%   See also GETDOF.

nNod=numel(unique(Nodenumbers(1:20)));

if nNod==15  % Degenerate prismatic element
  dof = dof_solid15(Nodenumbers([1 2 3 5 6 7 9 10 12 13 14 16 17 18 19]));
else
  dof=zeros(60,1);
  dof(1:3)  = Nodenumbers(1)  + [0.01 0.02 0.03]';
  dof(4:6)  = Nodenumbers(2)  + [0.01 0.02 0.03]';
  dof(7:9)  = Nodenumbers(3)  + [0.01 0.02 0.03]';
  dof(10:12)= Nodenumbers(4)  + [0.01 0.02 0.03]';
  dof(13:15)= Nodenumbers(5)  + [0.01 0.02 0.03]';
  dof(16:18)= Nodenumbers(6)  + [0.01 0.02 0.03]';
  dof(19:21)= Nodenumbers(7)  + [0.01 0.02 0.03]';
  dof(22:24)= Nodenumbers(8)  + [0.01 0.02 0.03]';
  dof(25:27)= Nodenumbers(9)  + [0.01 0.02 0.03]';
  dof(28:30)= Nodenumbers(10) + [0.01 0.02 0.03]';
  dof(31:33)= Nodenumbers(11) + [0.01 0.02 0.03]';
  dof(34:36)= Nodenumbers(12) + [0.01 0.02 0.03]';
  dof(37:39)= Nodenumbers(13) + [0.01 0.02 0.03]';
  dof(40:42)= Nodenumbers(14) + [0.01 0.02 0.03]';
  dof(43:45)= Nodenumbers(15) + [0.01 0.02 0.03]';
  dof(46:48)= Nodenumbers(16) + [0.01 0.02 0.03]';
  dof(49:51)= Nodenumbers(17) + [0.01 0.02 0.03]';
  dof(52:54)= Nodenumbers(18) + [0.01 0.02 0.03]';
  dof(55:57)= Nodenumbers(19) + [0.01 0.02 0.03]';
  dof(58:60)= Nodenumbers(20) + [0.01 0.02 0.03]';
end