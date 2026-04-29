function [X,Y,Z]=coord_rshell(Nodes,NodeNum)

%COORD_RSHELL  Coordinates of RSHELL element sides for plotting.
%
%   coord_rshell(Nodes,NodeNum)
%   returns the coordinates of RSHELL element sides for plotting.
%
%   Nodes       Node definitions [NodID x y z] (nNodes * 4)
%   NodeNum Node numbers [NodID1 NodID2 NodID3 NodID4] (nElem * 4)
%   X           X coordinates (4 * nElem)
%   Y           Y coordinates (4 * nElem)
%   Z           Z coordinates (4 * nElem)

% Mattias Schevenels
% April 2008

nElem=size(NodeNum,1);
X=zeros(4,nElem);
Y=zeros(4,nElem);
Z=zeros(4,nElem);

for iElem=1:nElem
  for iNode=1:4
    loc=find(Nodes(:,1)==NodeNum(iElem,iNode));
    if isempty(loc)
      error('Node %i is not defined.',NodeNum(iElem,iNode))
    elseif length(loc)>1
      error('Node %i is multiply defined.',NodeNum(iElem,iNode))
    else
      X(iNode,iElem)=Nodes(loc,2);
      Y(iNode,iElem)=Nodes(loc,3);
      Z(iNode,iElem)=Nodes(loc,4);
    end
  end
end
