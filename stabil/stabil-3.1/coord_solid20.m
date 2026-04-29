function [X,Y,Z]=coord_solid20(Nodes,Nodenumbers)

%COORD_SOLID20  Coordinates of SOLID20 element sides for plotting.
%
%   [X,Y,Z]=coord_solid20(Nodes,Nodenumbers)
%   returns the coordinates of SOLID20 element sides for plotting.
%
%   Nodes       Node definitions [NodID x y z] (nNodes * 4)
%   Nodenumbers Node numbers [NodID1 NodID2 NodID3 NodID4] (nElem * 4)
%   X           X coordinates (4 * nElem)
%   Y           Y coordinates (4 * nElem)
%   Z           Z coordinates (4 * nElem)

% David Dooms, Mattias Schevenels
% April 2008

nElem=size(Nodenumbers,1);
X=zeros(8,nElem);
Y=zeros(8,nElem);
Z=zeros(8,nElem);

for iElem=1:nElem
  for iNode=1:8
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

X=X([1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4],:);
Y=Y([1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4],:);
Z=Z([1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4],:);