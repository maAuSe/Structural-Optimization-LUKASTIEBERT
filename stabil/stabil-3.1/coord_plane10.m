function [X,Y,Z]=coord_plane15(Nodes,NodeNum)

%COORD_PLANE10  Coordinates of the plane elements for plotting.
%
%   [X,Y,Z] = coord_plane10(Nodes,NodeNum)
%   returns the coordinates of the plane10 elements for plotting.
%
%   Nodes      Node definitions        [NodID x y z] (nNodes * 4)
%   NodeNum    Node numbers            [NodID1 NodID2 ...] (nElem * 10)
%   X          X coordinates  (15 * nElem)    
%   Y          Y coordinates  (15 * nElem)
%   Z          Z coordinates  (15 * nElem)
%
%   See also COORD_TRUSS, PLOTELEM.

% Miche Jansen
% 2017


nNode=10;
nPoints=5;

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
  
  % Local coordinates of edges
  s = linspace(0,1,nPoints).';
  s_= s(end:-1:1);
  st=[s                zeros(nPoints,1) 
      s_               s
      zeros(nPoints,1) s_ ];
  
  N = zeros(3*(nPoints),nNode);
  for ind = 1:3*(nPoints)
      Ni = sh_t10(st(ind,1),st(ind,2));
      N(ind,:) = Ni; 
  end
  
  X = N*Xn;
  Y = N*Yn;
  Z = N*Zn;
else 
  X = [];
  Y = []; 
  Z = [];
end

end