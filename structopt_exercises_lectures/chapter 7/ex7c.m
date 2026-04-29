%EX7C
%
%   Shape optimization of a plane sheet for minimum compliance, using an
%   unstructured mesh.

% Mattias Schevenels
% January 2026

% MATERIAL, GEOMETRY, MESH, AND LOADING PARAMETERS
t = 0.01;        % Sheet thickness [m]
E = 200e9;       % Young's modulus [N/m^2]
nu = 0.3;        % Poisson's ratio
W = 0.85;        % Domain width [m]
H = 0.75;        % Initial domain height [m]
Vmax = 0.45*t;   % Maximum allowable volume [m^3]
P = -5e3;        % Point load [N]
elemSize = W/25; % Target maximum element side length [m]
p = 2;           % Degree of spline representing top boundary
n = 6;           % Number of control points along this spline

% DESIGN VARIABLES (y-coordinates of control points along top boundary, ordered left -> right)
x = repmat(H,n,1);
xmin = H/5;
xmax = H*2;

% ELEMENT TYPE, SECTION, AND MATERIAL PROPERTIES
Types = {1 'plane6'};
Sections = [1 t];
dSectionsdx = [];
Materials = [1 E nu];

% OPTIMIZATION LOOP
iter = 1;                     % Iteration counter
itermax = 1000;               % Maximum number of iterations
change = inf;                 % (Fictitious) initial design change
tol = 1e-4;                   % Convergence threshold
mmaparams = [];               % Internal MMA parameters
while true                   

  % CONTROL POINT COORDINATES (t = top / l = left / b = bottom / r = right boundary)
  Xt = linspace(0,W,n)';    Yt = x;                
  Xl = [0;0];               Yl = [0;x(1)];         
  Xb = [0;W];               Yb = [0;0];            
  Xr = [W;W];               Yr = [0;x(end)];       
  
  % CONTROL POINT COORDINATE SENSITIVITIES
  dXtdx = zeros(n,n);       dYtdx = eye(n);          
  dXldx = zeros(2,n);       dYldx = zeros(2,n);   dYldx(2,1)   = 1;
  dXbdx = zeros(2,n);       dYbdx = zeros(2,n);
  dXrdx = zeros(2,n);       dYrdx = zeros(2,n);   dYrdx(2,end) = 1;

  % CORNERS (used to estimate boundary lengths and number of nodes / sampling points per boundary)
  topleft = [0 x(1)];
  topright = [W x(end)];
  bottomleft = [0 0];
  bottomright = [W 0];

  % SPLINE SAMPLING POINTS (ordered in CCW direction, for Tracy's algorithm) 
  ut = linspace(1,0,ceil(norm(topleft-topright)/elemSize)*2+1);                         % Top: nodes ordered from right to left
  ul = linspace(1,0,ceil(norm(bottomleft-topleft)/elemSize)*2+1); ul = ul(2:end-1);     % Left: nodes ordered from top to bottom
  ub = linspace(0,1,ceil(norm(bottomright-bottomleft)/elemSize)*2+1);                   % Bottom: nodes ordered left to right
  ur = linspace(0,1,ceil(norm(topright-bottomright)/elemSize)*2+1); ur = ur(2:end-1);   % Right: nodes ordered bottom to top

  % SPLINE BASIS FUNCTIONS
  Bt = bspline(p,n,ut);
  Bl = bspline(1,2,ul);
  Bb = bspline(1,2,ub);
  Br = bspline(1,2,ur);

  % BOUNDARY NODE COORDINATES
  xt = Bt*Xt;               yt = Bt*Yt;              
  xb = Bb*Xb;               yb = Bb*Yb;             
  xl = Bl*Xl;               yl = Bl*Yl;             
  xr = Br*Xr;               yr = Br*Yr;   

  % BOUNDARY NODE COORDINATE SENSITIVITIES
  dxtdx = Bt*dXtdx;         dytdx = Bt*dYtdx;
  dxldx = Bl*dXldx;         dyldx = Bl*dYldx;
  dxbdx = Bb*dXbdx;         dybdx = Bb*dYbdx;
  dxrdx = Br*dXrdx;         dyrdx = Br*dYrdx;

  % BOUNDARY NODES
  xBoundaryNode = [xt; xl; xb; xr];    
  yBoundaryNode = [yt; yl; yb; yr];  
  nBoundaryNode = numel(xBoundaryNode);
  BoundaryNodes = [[1:nBoundaryNode]', xBoundaryNode, yBoundaryNode, zeros(nBoundaryNode,1)];

  % BOUNDARY NODE COORDINATE SENSITIVITIES
  dxBoundaryNodedx = [reshape(dxtdx,[],1,n);
                      reshape(dxldx,[],1,n);
                      reshape(dxbdx,[],1,n);
                      reshape(dxrdx,[],1,n)];
  dyBoundaryNodedx = [reshape(dytdx,[],1,n);
                      reshape(dyldx,[],1,n);
                      reshape(dybdx,[],1,n);
                      reshape(dyrdx,[],1,n)];
  dBoundaryNodesdx = zeros(nBoundaryNode,4,n);
  dBoundaryNodesdx(:,2,:) = dxBoundaryNodedx;
  dBoundaryNodesdx(:,3,:) = dyBoundaryNodedx;

  % INTERIOR NODES AND ELEMENTS
  [Nodes,Elements] = tracymesh6(BoundaryNodes,[],1,1,1);
  dNodesdx = zeros([size(Nodes),n]);
  dNodesdx(1:nBoundaryNode,:,:) = dBoundaryNodesdx;

  % BOUNDARY CONDITIONS AND LOADS (identified by coordinates as we use an unstructured mesh)
  DOF = getdof(Elements,Types);
  [~,jNode] = ismember(floor(DOF),Nodes(:,1));
  xDOF = Nodes(jNode,2);
  yDOF = Nodes(jNode,3);
  dDOF = round(100*(DOF-floor(DOF)));
  keepDOF = xDOF~=0;
  DOF = DOF(keepDOF);
  xDOF = xDOF(keepDOF);
  yDOF = yDOF(keepDOF);
  dDOF = dDOF(keepDOF);
  F = zeros(size(DOF));
  F((xDOF==W)&(yDOF==0)&(dDOF==2)) = P;

  % OBJECTIVE FUNCTION AND CONSTRAINT
  [K,~,dKdx] = asmkm(Nodes,Elements,Types,Sections,Materials,DOF,dNodesdx,dSectionsdx);
  U = K\F;
  C = F'*U;
  Ve = elemvolumes(Nodes,Elements,Types,Sections);
  V = sum(Ve);
  f0 = C;
  f = V/Vmax-1;

  % SENSITIVITIES
  dCdx = zeros(1,n);
  for k = 1:n
    dCdx(k) = -U'*dKdx{k}*U;
  end
  [~,dVedx] = elemvolumes(Nodes,Elements,Types,Sections,dNodesdx,dSectionsdx);
  dVdx = sum(dVedx);
  df0dx = dCdx;
  dfdx = dVdx/Vmax;

  % MMA UPDATE
  [xnew,y,z,lmult,mmaparams,subp,change,history] = mma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams);

  % PLOT CURRENT DESIGN
  plotelem(Nodes,Elements,Types,'Numbering','off','GCS','off');
  hold('on');
  plotnodes(Nodes,'Numbering','off','GCS','off');
  plot([Xt;Xl;Xb;Xr],[Yt;Yl;Yb;Yr],'ro','MarkerFaceColor','r','MarkerSize',4);
  hold('off');
  xlim([0 W]);
  ylim([0 H]);
  drawnow;

  % CHECK CONVERGENCE CRITERIA AND MOVE ON TO NEXT ITERATION
  if change>tol && (iter<itermax || isinf(change))
    iter = iter+1;
    x = xnew;
  else
    break;
  end
end
