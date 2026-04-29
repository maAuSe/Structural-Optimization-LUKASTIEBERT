function [Nodes,Elements] = tracymesh15(Nodes0,Elements0,Type,Section,Material,varargin)

%TRACYMESH15   Generate plane15 finite element mesh.
%
%   [Nodes,Elements] = TRACYMESH15(Nodes,Elements,Type,Section,Material,Front)
%   generates a plane 15-node finite element mesh by means of an advancing front
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
%               The number of front nodes must be a multiple of 4.
%               Default: Nodes(:,1).
%
%   OUTPUT ARGUMENTS
%
%   Nodes       Predefined + generated nodes (nNode * 4).
%   Elements    Element definitions (nElem * 19).

% Mattias Schevenels
% November 2023

% DEFAULT INPUT ARGUMENTS
if length(varargin)==0 || isstr(varargin{1})
  Front = Nodes0(:,1);
else
  Front = varargin{1};
end
plotmesh = any(strcmpi(varargin,'PlotMesh'));
if plotmesh sPlotmesh = 'PlotMesh'; else sPlotmesh=''; end

% CHECK IF FRONT IS SUITABLE FOR PLANE15 ELEMENTS
if rem(length(Front),4)~=0
  error('TRACYMESH15 expects the number of nodes on the front of the domain to mesh to be a multiple of 4.');
end

% GENERATE PLANE3 MESH
[Nodes,Elements1] = tracymesh3(Nodes0,[],Type,Section,Material,Front(1:4:end),sPlotmesh);
Elements1 = [Elements1 zeros(size(Elements1,1),12)];

% SUBDIVIDE PLANE3 ELEMENT BOUNDARIES
[int1,int2,int3] = deal(sparse([],[],[],max(Nodes(:,1)),max(Nodes(:,1)),3*size(Elements1,1)));
for k = 1:4:length(Front)
  n1 = Front(k);
  n2 = Front(rem(k+4,length(Front)));
  int1(n1,n2) = Front(k+1);
  int2(n1,n2) = Front(k+2);
  int3(n1,n2) = Front(k+3);
  int1(n2,n1) = Front(k+3);
  int2(n2,n1) = Front(k+2);
  int3(n2,n1) = Front(k+1);
end
lambda = [1 0 0;0 1 0;0 0 1;3/4 1/4 0;1/2 1/2 0;1/4 3/4 0;0 3/4 1/4;0 1/2 1/2;0 1/4 3/4;1/4 0 3/4;1/2 0 1/2;3/4 0 1/4;1/2 1/4 1/4;1/4 1/2 1/4;1/4 1/4 1/2];
lambda1 = lambda(4:6,:); % boundary 1
lambda2 = lambda(7:9,:); % boundary 2
lambda3 = lambda(10:12,:); % boundary 3
lambda4 = lambda(13:15,:); % inner

for k = 1:size(Elements1,1)
  n1 = Elements1(k,5);
  n2 = Elements1(k,6);
  n3 = Elements1(k,7);
  X = Nodes([n1 n2 n3],2:4);

  % boundary 1
  if int1(n1,n2)==0
    x = lambda1*X;
    newNodeNrs = max(Nodes(:,1))+[1;2;3];
    Nodes = [Nodes; newNodeNrs, x];
    int1(n1,n2) = newNodeNrs(1);
    int2(n1,n2) = newNodeNrs(2);
    int3(n1,n2) = newNodeNrs(3);
    int1(n2,n1) = newNodeNrs(3);
    int2(n2,n1) = newNodeNrs(2);
    int3(n2,n1) = newNodeNrs(1);
  end
  Elements1(k,8) = int1(n1,n2);
  Elements1(k,9) = int2(n1,n2);
  Elements1(k,10) = int3(n1,n2);

  % boundary 2
  if int1(n2,n3)==0
    x = lambda2*X;
    newNodeNrs = max(Nodes(:,1))+[1;2;3];
    Nodes = [Nodes; newNodeNrs, x];
    int1(n2,n3) = newNodeNrs(1);
    int2(n2,n3) = newNodeNrs(2);
    int3(n2,n3) = newNodeNrs(3);
    int1(n3,n2) = newNodeNrs(3);
    int2(n3,n2) = newNodeNrs(2);
    int3(n3,n2) = newNodeNrs(1);
  end
  Elements1(k,11) = int1(n2,n3);
  Elements1(k,12) = int2(n2,n3);
  Elements1(k,13) = int3(n2,n3);

  % boundary 3
  if int1(n3,n1)==0
    x = lambda3*X;
    newNodeNrs = max(Nodes(:,1))+[1;2;3];
    Nodes = [Nodes; newNodeNrs, x];
    int1(n3,n1) = newNodeNrs(1);
    int2(n3,n1) = newNodeNrs(2);
    int3(n3,n1) = newNodeNrs(3);
    int1(n1,n3) = newNodeNrs(3);
    int2(n1,n3) = newNodeNrs(2);
    int3(n1,n3) = newNodeNrs(1);
  end
  Elements1(k,14) = int1(n3,n1);
  Elements1(k,15) = int2(n3,n1);
  Elements1(k,16) = int3(n3,n1);

  % inner nodes
  x = lambda4*X;
  newNodeNrs = max(Nodes(:,1))+[1;2;3];
  Nodes = [Nodes; newNodeNrs, x];
  Elements1(k,17) = newNodeNrs(1);
  Elements1(k,18) = newNodeNrs(2);
  Elements1(k,19) = newNodeNrs(3);




end

% ADD NEW ELEMENTS TO PRE-EXISTING ELEMENTS
Elements = zeros(size(Elements0,1)+size(Elements1,1),max(size(Elements0,2),size(Elements1,2)));
Elements(1:size(Elements0,1),1:size(Elements0,2)) = Elements0;
Elements(size(Elements0,1)+[1:size(Elements1,1)],1:size(Elements1,2)) = Elements1;
