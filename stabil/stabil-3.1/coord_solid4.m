function [X,Y,Z]=coord_solid4(Nodes,Nodenumbers)

%COORD_SOLID4  Coordinates of SOLID4 element sides for plotting.
%
%   [X,Y,Z]=coord_rshell(Nodes,Nodenumbers)
%   returns the coordinates of RSHELL element sides for plotting.
%
%   Nodes       Node definitions [NodID x y z] (nNodes * 4)
%   Nodenumbers Node numbers [NodID1 NodID2 NodID3 NodID4] (nElem * 4)
%   X           X coordinates (4 * nElem)
%   Y           Y coordinates (4 * nElem)
%   Z           Z coordinates (4 * nElem)

nElem=size(Nodenumbers,1);
X=zeros(4,nElem);
Y=zeros(4,nElem);
Z=zeros(4,nElem);

for iElem=1:nElem
  for iNode=1:4
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

X=X([1 2 3 1 4 2 3 4],:);
Y=Y([1 2 3 1 4 2 3 4],:);
Z=Z([1 2 3 1 4 2 3 4],:);