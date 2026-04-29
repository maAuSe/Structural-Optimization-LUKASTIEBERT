function [pxyz,pind,pvalue]=patch_solid4(Nodes,NodeNum,Values)

%PATCH_SOLID4  Patch information of the solid4 elements for plotting.
%
%   [pxyz,pind,pvalue] = patch_solid20(Nodes,NodeNum,Values) 
%   returns matrices to plot patches of solid4 elements.
%
%   Nodes      Node definitions        [NodID x y z]
%   NodeNum    Node numbers       [NodID1 NodID2 NodID3 NodID4] (nElem * 10)
%   Values     Values assigned to nodes used for coloring    (nElem * 4)
%   pxyz       Coordinates of Nodes                          (4*nElem * 3)
%   pind       indices of Nodes                                (nElem * 4)
%   pvalue     Values arranged per Node                      (4*nElem * 1)
%
%   See also PATCH_SHELL20.

nElem=size(NodeNum,1);
NodeNum = NodeNum(:,1:4).'; % corner nodes only
if size(Values,2)< 4
  Values=Values(:,ones(1,4)).';    
else
  Values=Values(:,1:4).';
end

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
pvalue = Values(:);

block = [1 2 3;1 2 4;2 3 4;3 1 4]; % indices of faces of 1 block
pind = repmat(block,nElem,1)+ reshape(repmat(4*(0:(nElem-1)),12,1),3,4*nElem).';

% only draw outer (unique) faces ... 
blockt = block.';
[aa,aa,bb] = unique(sort(reshape(NodeNum(blockt(:),:),3,4*nElem).',2),'rows');
count = accumarray(bb(:),1);
pind = pind(count(bb)==1,:);
% reduce pxyz, pvalue?
[ind,tmp,pind] = unique(pind(:));
pind = reshape(pind,[],3);
pxyz = pxyz(ind,:);
pvalue = pvalue(ind,:);
end