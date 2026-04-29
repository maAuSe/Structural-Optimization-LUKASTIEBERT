function [Nodes,Elements] = tracymesh6(Nodes0,Elements0,Type,Section,Material,varargin)

%TRACYMESH6   Generate plane6 finite element mesh.
%
%   [Nodes,Elements] = TRACYMESH6(Nodes,Elements,Type,Section,Material,Front)
%   generates a plane 6-node finite element mesh by means of an advancing front
%   algorithm (Tracy, 1977) followed by a regularization step based on the force
%   density method (Schek, 1974).  The result is appended to the Nodes and the
%   Elements matrices.
%
%   INPUT ARGUMENTS
%
%   Nodes       Pre-existing nodes, including all nodes on the boundary of the
%               area to mesh, defined following Stabil's convention (nNode0 * 4).
%   Elements    Pre-existing elements, defined following Stabil's convention (nElem0 * x).
%   Type        Element type identifier used to populate the 2nd column of
%               the output Elements matrix (1 * 1).
%   Section     Section identifier (1 * 1).
%   Material    Material identifier (1 * 1).
%   Front       Identifiers of the nodes on the boundary of the area to mesh.
%               This area must be located in the XY-plane.
%               Nodes on the boundary must be sorted in counterclockwise order.
%               The number of front nodes must be a multiple of 2.
%               Default: Nodes(:,1).
%
%   OUTPUT ARGUMENTS
%
%   Nodes       Predefined + generated nodes (nNode * 4).
%   Elements    Element definitions (nElem * 10).

% Mattias Schevenels
% January 2024

% DEFAULT INPUT ARGUMENTS
if length(varargin)==0 || isstr(varargin{1})
  Front = Nodes0(:,1);
else
  Front = varargin{1};
end
plotmesh = any(strcmpi(varargin,'PlotMesh'));
if plotmesh sPlotmesh = 'PlotMesh'; else sPlotmesh=''; end

% CHECK IF FRONT IS SUITABLE FOR PLANE6 ELEMENTS
if rem(length(Front),2)~=0
  error('TRACYMESH6 expects the number of nodes on the front of the domain to mesh to be even.');
end

% GENERATE PLANE3 MESH
[Nodes,Elements1] = tracymesh3(Nodes0,[],Type,Section,Material,Front(1:2:end),sPlotmesh);
Elements1 = [Elements1 zeros(size(Elements1,1),3)];

% SUBDIVIDE PLANE3 ELEMENT BOUNDARIES
int = sparse([],[],[],max(Nodes(:,1)),max(Nodes(:,1)),size(Elements1,1));
for k = 1:2:length(Front)
  n1 = Front(k);
  n2 = Front(rem(k+2,length(Front)));
  int(n1,n2) = Front(k+1);
  int(n2,n1) = Front(k+1);
end
lambda = [1 0 0;0 1 0;0 0 1;1/2 1/2 0;0 1/2 1/2;1/2 0 1/2];
lambda1 = lambda(4,:); % boundary 1
lambda2 = lambda(5,:); % boundary 2
lambda3 = lambda(6,:); % boundary 3

for k = 1:size(Elements1,1)
  n1 = Elements1(k,5);
  n2 = Elements1(k,6);
  n3 = Elements1(k,7);
  X = Nodes([n1 n2 n3],2:4);

  % boundary 1
  if int(n1,n2)==0
    x = lambda1*X;
    newNodeNr = max(Nodes(:,1))+1;
    Nodes = [Nodes; newNodeNr, x];
    int(n1,n2) = newNodeNr;
    int(n2,n1) = newNodeNr;
  end
  Elements1(k,8) = int(n1,n2);

  % boundary 2
  if int(n2,n3)==0
    x = lambda2*X;
    newNodeNr = max(Nodes(:,1))+1;
    Nodes = [Nodes; newNodeNr, x];
    int(n2,n3) = newNodeNr;
    int(n3,n2) = newNodeNr;
  end
  Elements1(k,9) = int(n2,n3);

  % boundary 3
  if int(n3,n1)==0
    x = lambda3*X;
    newNodeNr = max(Nodes(:,1))+1;
    Nodes = [Nodes; newNodeNr, x];
    int(n3,n1) = newNodeNr;
    int(n1,n3) = newNodeNr;
  end
  Elements1(k,10) = int(n3,n1);

end

% ADD NEW ELEMENTS TO PRE-EXISTING ELEMENTS
Elements = zeros(size(Elements0,1)+size(Elements1,1),max(size(Elements0,2),size(Elements1,2)));
Elements(1:size(Elements0,1),1:size(Elements0,2)) = Elements0;
Elements(size(Elements0,1)+[1:size(Elements1,1)],1:size(Elements1,2)) = Elements1;
