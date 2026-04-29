function [pxyz,pind,pvalue]=patch_mass(Nodes,NodeNum,Values)
%PATCH_MASS  Patch information of the mass elements for plotting.
%
%   [pxyz,pind,pvalue] = patch_mass(Nodes,NodeNum,Values) returns matrices
%   to plot patches of mass elements.
%
%   Nodes      Node definitions [NodID x y z].
%   NodeNum    Node numbers [NodID1 NodID2 NodID3] (nElem).
%   Values     Values assigned to nodes used for coloring (nElem).
%   pxyz       Coordinates of Nodes  (nElem * 3).
%   pind       Indices of Nodes  (nElem).
%   pvalue     Values arranged per Node (nElem).

nElem=size(NodeNum,1);

pxyz = zeros(1*nElem,3);


for iElem=1:nElem
    for ind = 1:1  
    loc=find(Nodes(:,1)==NodeNum(iElem,ind));
    if isempty(loc)
        error('Node %i is not defined.',NodeNum(iElem,1))
    elseif length(loc)>1
        error('Node %i is multiply defined.',NodeNum(iElem,1))
    else
        pxyz(1*(iElem-1)+ind,1)=Nodes(loc,2);
        pxyz(1*(iElem-1)+ind,2)=Nodes(loc,3);
        pxyz(1*(iElem-1)+ind,3)=Nodes(loc,4);
    end
    end
   
end

pind = reshape((1:2*nElem),2,nElem);
pind = pind.';

pvalue = nan(size(pxyz,1),1);