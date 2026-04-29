function [X,Y,Z]=coord_shell8(Nodes,NodeNum)

%COORD_SHELL8  Coordinates of the shell8 elements for plotting.
%
%   [X,Y,Z] = coord_shell8(Nodes,NodeNum)
%   returns the coordinates of the shell8 elements for plotting.
%
%   Nodes      Node definitions        [NodID x y z] (nNodes * 4)
%   NodeNum    Node numbers            [NodID1 NodID2 NodID3 NodID4] (nElem * 8)
%   X          X coordinates  (20 * nElem)    
%   Y          Y coordinates  (20 * nElem)
%   Z          Z coordinates  (20 * nElem)
%
%   See also COORD_TRUSS, PLOTELEM.

% Miche Jansen
% 2009

nNode = 8; % Number of nodes per element
nPoints = 5; % aanpassen voor betere weergave (indien element gekromd is)

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

N = zeros(4*(nPoints),8);
s = linspace(-1,1,nPoints).';
st = [s -ones(nPoints,1);
      ones(nPoints,1) s;
      -s ones(nPoints,1);
      -ones(nPoints,1) -s];


for ind = 1:4*(nPoints)
    Ni = sh_qs8(st(ind,1),st(ind,2));
    N(ind,:) = Ni; 
end

X = N*Xn;
Y = N*Yn;
Z = N*Zn;

else X = [];Y = []; Z = [];
end

end