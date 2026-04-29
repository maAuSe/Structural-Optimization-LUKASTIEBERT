function [X,Y,Z]=coord_shell6(Nodes,NodeNum)

%COORD_SHELL6  Coordinates of the shell6 elements for plotting.
%
%   [X,Y,Z] = coord_shell6(Nodes,NodeNum)
%   returns the coordinates of the shell6 elements for plotting.
%
%   Nodes      Node definitions        [NodID x y z] (nNodes * 4)
%   NodeNum    Node numbers            [NodID1 NodID2 NodID3 NodID4] (nElem * 6)
%   X          X coordinates  (15 * nElem)    
%   Y          Y coordinates  (15 * nElem)
%   Z          Z coordinates  (15 * nElem)
%
%   See also COORD_TRUSS, PLOTELEM.

% Miche Jansen
% 2009

nNode=6;   % Number of nodes per element
nPoints=5; % aanpassen voor betere weergave (indien element gekromd is)

% check nodes
if length(unique(Nodes(:,1)))~=length(Nodes(:,1))
    for iNode = 1:length(Nodes(:,1))
        loc = find(abs(Nodes(iNode,1)-Nodes(:,1))<1e-4);
        if length(loc)>1
        error('Node %i is multiply defined.',Nodes(iNode,1))
        end
    end
end

if ~isempty(NodeNum)

Xn = Nodes(:,2);
Yn = Nodes(:,3);
Zn = Nodes(:,4);

NodeNum = NodeNum(:,1:nNode);
[exis,loc] = ismember(round(NodeNum.'),round(Nodes(:,1)));

if ~all(exis(:))
   unknowns = NodeNum(~exis.');
   error('Node %i is not defined.',unknowns(1));
end

Xn = Xn(loc);
Yn = Yn(loc);
Zn = Zn(loc);

N = zeros(3*(nPoints),6);
s = linspace(0,1,nPoints).';
s_= s(end:-1:1);
st=[s                zeros(nPoints,1) 
    s_               s
    zeros(nPoints,1) s_ ];

for ind = 1:3*(nPoints)
    Ni = sh_t6(st(ind,1),st(ind,2));
    N(ind,:) = Ni; 
end

X = N*Xn;
Y = N*Yn;
Z = N*Zn;

else X = [];Y = []; Z = [];
end

end