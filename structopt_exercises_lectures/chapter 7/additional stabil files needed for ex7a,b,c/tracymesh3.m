function [Nodes,Elements] = tracymesh3(Nodes0,Elements0,Type,Section,Material,varargin)

%TRACYMESH3   Generate plane3 finite element mesh.
%
%   [Nodes,Elements] = TRACYMESH3(Nodes,Elements,Type,Section,Material,Front)
%   generates a plane 3-node finite element mesh by means of an advancing front
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
%               Default: Nodes(:,1).
%
%   OUTPUT ARGUMENTS
%
%   Nodes       Predefined + generated nodes (nNode * 4).
%   Elements    Element definitions (nElem * 7).

% Mattias Schevenels
% November 2023

% DEFAULT INPUT ARGUMENTS
if length(varargin)==0 || isstr(varargin{1})
  Front = Nodes0(:,1);
else
  Front = varargin{1};
end
plotmesh = any(strcmpi(varargin,'PlotMesh'));
if isempty(Elements0), Elements0 = zeros(0,7); end

% CHECK INPUT ARGUMENTS
if any(Nodes0(Front,4)~=0), error('Function TRACYMESH3 only works for 2D meshes in the XY-plane.'); end

% TRACY'S ADVANCING FRONT ALGORITHM
Nodes = Nodes0;
Elements1 = zeros(0,7);

k = Front;                                       % front node indices
x = Nodes(k,2);                                  % front node x-coordinates
y = Nodes(k,3);                                  % front node y-coordinates

dx1 = x-x([end,1:end-1]); dx2 = x([2:end,1])-x;  % difference in x-coordinates wrt previous/next node on front
dy1 = y-y([end,1:end-1]); dy2 = y([2:end,1])-y;  % difference in y-coordinates wrt previous/next node on front
theta = sign(dx1.*dy2-dx2.*dy1).*real(acos((dx1.*dx2+dy1.*dy2)./sqrt((dx1.^2+dy1.^2).*(dx2.^2+dy2.^2)))); % change of angle between two adjacent parts of the front

d = (sqrt(dx1.^2+dy1.^2)+sqrt(dx2.^2+dy2.^2))/2; % average distance wrt adjacent nodes on front
elemsize = scatteredInterpolant(Nodes0(k,2),Nodes0(k,3),d); % function to determine target element side lengths throughout mesh domain, obtained as an interpolation of the element side lengths on the boundary

% ESTIMATE UPPER BOUND FOR NUMBER OF TRIANGLES
perimeter = sum(sqrt(dx1.^2+dy1.^2));     % upper bound of perimeter
R = perimeter/pi/2;                       % radius if domain would be circular
areaD = pi*R^2;                           % upper bound of domain area
areaT = sqrt(3)/4*min(d)^2;               % area of smallest triangle (if equilateral)
maxNumTriangles = 5*areaD/areaT;          % maximum expected number of triangles (including safety factor)

% LOOP UNTIL FRONT REPRESENTS A SINGLE TRIANGLE
while length(k)>3

  % CREATE A NEW ELEMENT IN ALL CORNER ANGLES LESS THAN 90 DEGREES (STARTING WITH THE SHARPEST CORNER ANGLE)
  while length(k)>3 && any(theta>pi/2)
    [~,i1] = max(theta);
    i2 = i1+1; if i2>length(k), i2 = 1; end
    i0 = i1-1; if i0<1, i0 = length(k); end
    k0 = k(i0); k1 = k(i1); k2 = k(i2);

    % ADD NEW ELEMENT
    Elements1 = [Elements1; max([0;Elements0(:,1);Elements1(:,1)])+1, Type, Section, Material, k1 k2 k0];

    % UPDATE FRONT DEFINITION
    k = k([1:i1-1,i1+1:end]);
    x = x([1:i1-1,i1+1:end]);
    y = y([1:i1-1,i1+1:end]);
    dx1 = x-x([end,1:end-1]); dx2 = x([2:end,1])-x;
    dy1 = y-y([end,1:end-1]); dy2 = y([2:end,1])-y;
    theta = sign(dx1.*dy2-dx2.*dy1).*real(acos((dx1.*dx2+dy1.*dy2)./sqrt((dx1.^2+dy1.^2).*(dx2.^2+dy2.^2))));

    % CHECK NUMBER OF TRIANGLES AND THROW ERROR IF IT BECOMES SUSPICIOUSLY LARGE
    if size(Elements1,1)>maxNumTriangles
      figure;
      plotelem(Nodes,Elements1(2:end,:),{1 'plane3'},'Numbering','off');
      error('Maximum expected number of triangles exceeded; probably TRACYMESH3 ran into problems - check figure.');
    end

    % PLOT MESH
    if plotmesh
      plotelem(Nodes,Elements1(2:end,:),{1 'plane3'},'Numbering','off');
      drawnow;
    end
  end

  % DEFINE A NEW NODE AND CREATE TWO NEW ELEMENTS A CORNER ANGLE LESS THAN 180 DEGREES (THE SMALLEST CORNER ANGLE)
  if length(k)>3 && any(theta>0)
    [~,i1] = max(theta);
    i2 = i1+1; if i2>length(k), i2 = 1; end
    i0 = i1-1; if i0<1, i0 = length(k); end
    k0 = k(i0); k1 = k(i1); k2 = k(i2);
    x0 = x(i0); x1 = x(i1); x2 = x(i2);
    y0 = y(i0); y1 = y(i1); y2 = y(i2);
    x3 = (x0+x2)/2 + (y0-y2)/5;
    y3 = (y0+y2)/2 + (x2-x0)/5;
    k3 = max(Nodes(:,1))+1;

    % MODIFICATION WRT TRACY'S ALGORITHM: CHOOSE NEW NODE TAKING INTO ACCOUNT TARGET ELEMENT SIZE AT PRESENT LOCATION AND DISTANCE TO OTHER SIDE OF FRONT (AVOID OVERLAPS)
    d1 = sqrt((x0-x1)^2+(y0-y1)^2);
    d2 = sqrt((x2-x1)^2+(y2-y1)^2);
    d3 = sqrt((x3-x1)^2+(y3-y1)^2);
    iall = 1:length(k);
    iother = setdiff(iall,[i0,i1,i2]);
    dmax = 0.7*min(sqrt((x1-x(iother)).^2+(y1-y(iother)).^2));
    dtarget = min(dmax,elemsize(x3,y3));
    x3 = x1+(x3-x1)/d3*dtarget;
    y3 = y1+(y3-y1)/d3*dtarget;

    % ADD NEW NODE AND NEW ELEMENTS
    Nodes = [Nodes; k3, x3, y3, 0];
    Elements1 = [Elements1; max([0;Elements0(:,1);Elements1(:,1)])+1, Type, Section, Material, k1 k3 k0;
                            max([0;Elements0(:,1);Elements1(:,1)])+2, Type, Section, Material, k1 k2 k3];

    % UPDATE FRONT DEFINITION
    x(i1) = x3;
    y(i1) = y3;
    k(i1) = k3;
    dx1 = x-x([end,1:end-1]); dx2 = x([2:end,1])-x;
    dy1 = y-y([end,1:end-1]); dy2 = y([2:end,1])-y;
    theta = sign(dx1.*dy2-dx2.*dy1).*real(acos((dx1.*dx2+dy1.*dy2)./sqrt((dx1.^2+dy1.^2).*(dx2.^2+dy2.^2))));

    % PLOT MESH
    if plotmesh
      plotelem(Nodes,Elements1,{1 'plane3'},'Numbering','off');
      drawnow;
    end
  end
end

% ADD ELEMENT CORRESPONDING TO FINAL REMAINING TRIANGLE
Elements1 = [Elements1; max([0;Elements0(:,1);Elements1(:,1)])+1, Type, Section, Material, k'];

% ADD NEW ELEMENTS TO PRE-EXISTING ELEMENTS
Elements = zeros(size(Elements0,1)+size(Elements1,1),max(size(Elements0,2),size(Elements1,2)));
Elements(1:size(Elements0,1),1:size(Elements0,2)) = Elements0;
Elements(size(Elements0,1)+[1:size(Elements1,1)],1:size(Elements1,2)) = Elements1;

% REGULARIZE MESH USING FDM
n1n2 = unique([min(Elements1(:,[5,6]),[],2), max(Elements1(:,[5,6]),[],2);
               min(Elements1(:,[6,7]),[],2), max(Elements1(:,[6,7]),[],2);
               min(Elements1(:,[7,5]),[],2), max(Elements1(:,[7,5]),[],2)],'rows');
n1 = n1n2(:,1);
n2 = n1n2(:,2);
elem = [1:length(n1)];
nElem = length(elem);
%nNode = max([n1;n2]);
nNode = max(Nodes(:,1));
Cs = sparse(elem,n1,-1,nElem,nNode)+sparse(elem,n2,1,nElem,nNode);
fixednodes = Nodes0(:,1);
xf = Nodes0(:,2);
yf = Nodes0(:,3);
allnodes = [1:nNode]';
varnodes = setdiff(allnodes,fixednodes);
C = Cs(:,varnodes);
Cf = Cs(:,fixednodes);
D = C'*C;
Df = C'*Cf;
x = D\(-Df*xf);
y = D\(-Df*yf);
[xs,ys] = deal(zeros(nNode,1));
xs([varnodes;fixednodes]) = [x;xf];
ys([varnodes;fixednodes]) = [y;yf];
Nodes(:,2:3) = [xs ys];

% PLOT MESH
if plotmesh
  plotelem(Nodes,Elements1,{1 'plane3'},'Numbering','off');
  drawnow;
end
