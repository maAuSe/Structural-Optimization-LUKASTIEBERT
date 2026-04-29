function [X,Y,Z]=coord_plane3(Nodes,Nodenumbers)

%COORD_PLANE3  Coordinates of PLANE3 element sides for plotting.
%
%   [X,Y,Z]=coord_plane3(Nodes,Nodenumbers)
%   returns the coordinates of PLANE3 element sides for plotting.
%
%   Nodes       Node definitions [NodID x y z] (nNodes * 3)
%   Nodenumbers Node numbers [NodID1 NodID2 NodID3] (nElem * 3)
%   X           X coordinates (3 * nElem)
%   Y           Y coordinates (3 * nElem)
%   Z           Z coordinates (3 * nElem)

nElem=size(Nodenumbers,1);
X=zeros(3,nElem);
Y=zeros(3,nElem);
Z=zeros(3,nElem);

for iElem=1:nElem
  for iNode=1:3
    loc=find(Nodes(:,1)==Nodenumbers(iElem,iNode));
    if isempty(loc)
      error('Node %i is not defined.',Nodenumbers(iElem,iNode))
    elseif length(loc)>1
      error('Node %i is multiply defined.',Nodenumbers(iElem,iNode))
    else
      X(iNode,iElem)=Nodes(loc,2);
      Y(iNode,iElem)=Nodes(loc,3);
      Z(iNode,iElem)=Nodes(loc,4);
    end
  end
end
