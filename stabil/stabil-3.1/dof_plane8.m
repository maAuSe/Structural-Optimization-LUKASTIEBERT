function dof = dof_plane8(Nodenumbers)

dof=zeros(16,1);
dof(1:2)   = Nodenumbers(1)+[0.01;0.02];
dof(3:4)   = Nodenumbers(2)+[0.01;0.02];
dof(5:6)   = Nodenumbers(3)+[0.01;0.02];
dof(7:8)   = Nodenumbers(4)+[0.01;0.02];
dof(9:10)  = Nodenumbers(5)+[0.01;0.02];
dof(11:12) = Nodenumbers(6)+[0.01;0.02];
dof(13:14) = Nodenumbers(7)+[0.01;0.02];
dof(15:16) = Nodenumbers(8)+[0.01;0.02];