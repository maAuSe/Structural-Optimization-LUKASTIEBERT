function dof = dof_rshell(NodeNum)

% Mattias Schevenels
% April 2008

dof=zeros(24,1);
dof(1:6)=NodeNum(1)+[0.01:0.01:0.06]';
dof(7:12)=NodeNum(2)+[0.01:0.01:0.06]';
dof(13:18)=NodeNum(3)+[0.01:0.01:0.06]';
dof(19:24)=NodeNum(4)+[0.01:0.01:0.06]';
