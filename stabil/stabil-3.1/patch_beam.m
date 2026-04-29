function [pxyz,pind,pvalue]=patch_beam(Nodes,NodeNum,Values)
%PATCH_BEAM  Patch information of the beam elements for plotting.
%
%   [pxyz,pind,pvalue] = patch_beam(Nodes,NodeNum,Values) returns matrices
%   to plot patches of beam elements.
%
%   Nodes      Node definitions [NodID x y z].
%   NodeNum    Node numbers [NodID1 NodID2 NodID3] (nElem * 3).
%   Values     Values assigned to nodes used for coloring (nElem * 3).
%   pxyz       Coordinates of Nodes  (3*nElem * 3).
%   pind       Indices of Nodes  (nElem * 3).
%   pvalue     Values arranged per Node (3*nElem * 1).

nElem=size(NodeNum,1);
NodeNum = NodeNum(:,1:2).';
pxyz = zeros(2*nElem,3);


% check nodes
if length(unique(Nodes(:,1)))~=length(Nodes(:,1))
    for iNode = 1:length(Nodes(:,1))
        loc = find(abs(Nodes(iNode,1)-Nodes(:,1))<1e-4);
        if length(loc)>1
        error('Node %i is multiply defined.',Nodes(iNode,1))
        end
    end
end

[exis,loc] = ismember(round(NodeNum),round(Nodes(:,1)));

if ~all(exis(:))
   unknowns = NodeNum(~exis);
   error('Node %i is not defined.',unknowns(1));
end

pxyz = Nodes(loc(:),2:4);
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