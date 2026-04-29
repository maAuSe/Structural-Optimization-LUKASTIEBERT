function [X,Y,Z]=coord_truss(Nodes,NodeNum)

%COORD_TRUSS   Coordinates of the truss elements for plotting.
%
%   [X,Y,Z]=coord_truss(Nodes,NodeNum)
%   returns the coordinates of the truss elements for plotting.
%
%   Nodes      Node definitions        [NodID x y z] (nNodes * 4)
%   NodeNum Node numbers          [NodID1 NodID2] (nElem * 2)
%   X          X coordinates  (2 * nElem)    
%   Y          Y coordinates  (2 * nElem)
%   Z          Z coordinates  (2 * nElem)
%
%   See also COORD_BEAM, PLOTELEM.


% David Dooms
% March 2008

nElem=size(NodeNum,1);
X=zeros(2,nElem);
Y=zeros(2,nElem);
Z=zeros(2,nElem);

for iElem=1:nElem
    loc1=find(Nodes(:,1)==NodeNum(iElem,1));
    loc2=find(Nodes(:,1)==NodeNum(iElem,2));
    if isempty(loc1)
        error('Node %i is not defined.',NodeNum(iElem,1))
    elseif length(loc1)>1
        error('Node %i is multiply defined.',NodeNum(iElem,1))
    elseif isempty(loc2)
        error('Node %i is not defined.',NodeNum(iElem,2))
    elseif length(loc2)>1
        error('Node %i is multiply defined.',NodeNum(iElem,2))
    else
        X(1,iElem)=Nodes(loc1,2);
        X(2,iElem)=Nodes(loc2,2);
        Y(1,iElem)=Nodes(loc1,3);
        Y(2,iElem)=Nodes(loc2,3);
        Z(1,iElem)=Nodes(loc1,4);
        Z(2,iElem)=Nodes(loc2,4);
    end
end   
