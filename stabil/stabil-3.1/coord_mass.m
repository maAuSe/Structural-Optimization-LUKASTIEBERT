function [X,Y,Z]=coord_mass(Nodes,NodeNum)

%COORD_MASS   Coordinates of the mass element for plotting.
%
%   [X,Y,Z]=coord_mass(Nodes,NodeNum)
%   returns the coordinates of the mass element for plotting.
%
%   Nodes      Node definitions        [NodID x y z] (nNode * 4)
%   NodeNum Node numbers               [NodID1] (nElem * 1)
%   X          X coordinates  (1 * nElem)    
%   Y          Y coordinates  (1 * nElem)
%   Z          Z coordinates  (1 * nElem)
%
%   See also COORD_BEAM, PLOTELEM.

nElem=size(NodeNum,1);
X=zeros(1,nElem);
Y=zeros(1,nElem);
Z=zeros(1,nElem);

for iElem=1:nElem
  loc1=find(Nodes(:,1)==NodeNum(iElem,1));
  if isempty(loc1)
    error('Node %i is not defined.',NodeNum(iElem,1))
  elseif length(loc1)>1
    error('Node %i is multiply defined.',NodeNum(iElem,1))
  else
    X(1,iElem)=Nodes(loc1,2);
    Y(1,iElem)=Nodes(loc1,3);
    Z(1,iElem)=Nodes(loc1,4);
  end
end   
