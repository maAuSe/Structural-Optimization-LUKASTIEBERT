function Nodesel=selectnode(Nodes,xmin,ymin,zmin,xmax,ymax,zmax)

%SELECTNODE   Select nodes by location.
%
%   Nodesel=selectnode(Nodes,x,y,z)
%   Nodesel=selectnode(Nodes,xmin,ymin,zmin,xmax,ymax,zmax)
%   selects nodes by location.
%
%   Nodes      Node definitions          [NodID x y z]
%   Nodesel    Node definitions of the selected nodes

% David Dooms
% March 2008

if nargin < 5 
    toler=1e-10*max(max(abs(Nodes(:,2:4))));   
    xmax=xmin+toler;
    xmin=xmin-toler;
    ymax=ymin+toler;
    ymin=ymin-toler;
    zmax=zmin+toler;
    zmin=zmin-toler;
end
     
Nodesel=Nodes((Nodes(:,2)>=xmin) & (Nodes(:,2)<=xmax) & ...
              (Nodes(:,3)>=ymin) & (Nodes(:,3)<=ymax) & ... 
              (Nodes(:,4)>=zmin) & (Nodes(:,4)<=zmax),:);
end
