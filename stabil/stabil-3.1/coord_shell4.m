function [X,Y,Z]=coord_shell4(Nodes,NodeNum)

%COORD_SHELL4  Coordinates of the shell elements for plotting.
%
%   [X,Y,Z] = coord_shell4(Nodes,NodeNum)
%   returns the coordinates of the shell4 elements for plotting.
%
%   Nodes      Node definitions        [NodID x y z] (nNodes * 4)
%   NodeNum    Node numbers            [NodID1 NodID2 NodID3 NodID4] (nElem * 4)
%   X          X coordinates  (4 * nElem)    
%   Y          Y coordinates  (4 * nElem)
%   Z          Z coordinates  (4 * nElem)
%
%   See also COORD_TRUSS, PLOTELEM.

% Miche Jansen
% 2009

% check nodes
if length(unique(Nodes(:,1)))~=length(Nodes(:,1))
    for iNode = 1:length(Nodes(:,1))
        loc = find(abs(Nodes(iNode,1)-Nodes(:,1))<1e-4);
        if length(loc)>1
        error('Node %i is multiply defined.',Nodes(iNode,1))
        end
    end
end

nNode = 4; % Number of nodes per element

X = Nodes(:,2);
Y = Nodes(:,3);
Z = Nodes(:,4);

NodeNum = NodeNum(:,1:nNode);
[exis,loc] = ismember(round(NodeNum.'),round(Nodes(:,1)));

if ~all(exis(:))
   unknowns = NodeNum(~exis.');
   error('Node %i is not defined.',unknowns(1));
end

X = X(loc);
Y = Y(loc);
Z = Z(loc);