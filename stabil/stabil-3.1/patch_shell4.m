function [pxyz,pind,pvalue]=patch_shell4(Nodes,NodeNum,Values)

%PATCH_SHELL4  Patch information of the shell4 elements for plotting.
%
%   [pxyz,pind,pvalue] = patch_shell4(Nodes,NodeNum,Values) returns matrices
%   to plot patches of shell4 elements.
%
%   Nodes      Node definitions        [NodID x y z]
%   NodeNum    Node numbers       [NodID1 NodID2 NodID3 NodID4] (nElem * 4)
%   Values     Values assigned to nodes used for coloring    (nElem * 4)
%   pxyz       Coordinates of Nodes                          (4*nElem * 3)
%   pind       indices of Nodes                                (nElem * 4)
%   pvalue     Values arranged per Node                      (4*nElem * 1)
%
%   See also PLOTSTRESSCONTOURF, PLOTSHELLFCONTOURF.

% Miche Jansen
% 2010

nElem=size(NodeNum,1);
NodeNum = NodeNum(:,1:4).';
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

pind = reshape((1:4*nElem),4,nElem);
pind = pind.';

end
