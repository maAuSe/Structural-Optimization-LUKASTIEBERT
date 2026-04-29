function [pxyz,pind,pvalue]=patch_truss(Nodes,NodeNum,Values)
%PATCH_TRUSS  Patch information of the truss elements for plotting.
%
%   [pxyz,pind,pvalue] = patch_truss(Nodes,NodeNum,Values) returns matrices
%   to plot patches of truss elements.
%
%   Nodes      Node definitions [NodID x y z].
%   NodeNum    Node numbers [NodID1 NodID2 NodID3] (nElem * 2).
%   Values     Values assigned to nodes used for coloring (nElem * 2).
%   pxyz       Coordinates of Nodes  (2*nElem * 3).
%   pind       Indices of Nodes  (nElem * 2).
%   pvalue     Values arranged per Node (2*nElem * 1).

nElem=size(NodeNum,1);

pxyz = zeros(2*nElem,3);


for iElem=1:nElem
    for ind = 1:2  
    loc=find(Nodes(:,1)==NodeNum(iElem,ind));
    if isempty(loc)
        error('Node %i is not defined.',NodeNum(iElem,1))
    elseif length(loc)>1
        error('Node %i is multiply defined.',NodeNum(iElem,1))
    else
        pxyz(2*(iElem-1)+ind,1)=Nodes(loc,2);
        pxyz(2*(iElem-1)+ind,2)=Nodes(loc,3);
        pxyz(2*(iElem-1)+ind,3)=Nodes(loc,4);
    end
    end
   
end

pind = reshape((1:2*nElem),2,nElem);
pind = pind.';

if isempty(Values)
    pvalue = nan(size(pxyz,1),1);
else
    if size(Values,2)< 2
    Values=Values(:,ones(1,2)).';    
    else
    Values=Values(:,1:2).';
    end
    pvalue = Values(:);
end