%EX7B
%
%   Shape optimization of a plane sheet for minimum compliance, using a
%   Coons surface mesh.

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
nelx = 30;       % Number of elements in x-direction
nely = 10;       % Number of elements in y-direction
nnx = 2*nelx+1;  % Number of nodes in x-direction
nny = 2*nely+1;  % Number of nodes in y-direction
p = 2;           % Degree of spline representing top boundary
n = 6;           % Number of control points along this spline

% DESIGN VARIABLES (y-coordinates of control points along top boundary, ordered left -> right)
x = repmat(H,n,1);
xmin = H/5;
xmax = H*2;

% SPLINE SAMPLING POINTS (top+bottom / left+right boundary)
u = linspace(0,1,nnx);
v = linspace(0,1,nny);

% SPLINE BASIS FUNCTIONS (top / left / bottom / right boundary)
Bt = bspline(p,n,u);          
Bl = bspline(1,2,v);          
Bb = bspline(1,2,u);          
Br = bspline(1,2,v);          

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

  % BOUNDARY NODES (node IDs assigned inside COONSMESH, first column is placeholder)
  NodesT = [zeros(nnx,1), xt, yt, zeros(nnx,1)];
  NodesL = [zeros(nny,1), xl, yl, zeros(nny,1)];
  NodesB = [zeros(nnx,1), xb, yb, zeros(nnx,1)];
  NodesR = [zeros(nny,1), xr, yr, zeros(nny,1)];

  % BOUNDARY NODE SENSITIVITIES
  dNodesTdx = zeros(nnx,4,n);
  dNodesLdx = zeros(nny,4,n);
  dNodesBdx = zeros(nnx,4,n);
  dNodesRdx = zeros(nny,4,n);

  dNodesTdx(:,2,:) = reshape(dxtdx, nnx, 1, n);
  dNodesTdx(:,3,:) = reshape(dytdx, nnx, 1, n);

  dNodesLdx(:,2,:) = reshape(dxldx, nny, 1, n);
  dNodesLdx(:,3,:) = reshape(dyldx, nny, 1, n);

  dNodesRdx(:,2,:) = reshape(dxrdx, nny, 1, n);
  dNodesRdx(:,3,:) = reshape(dyrdx, nny, 1, n);

  % COONS SURFACE
  [Nodes,dNodesdx,NodesT,NodesL,NodesB,NodesR] = coonsmesh(NodesT,NodesL,NodesB,NodesR,dNodesTdx,dNodesLdx,dNodesBdx,dNodesRdx);

  % ELEMENTS
  Elements = [1 1 1 1 1 3 3+2*nnx  2 3+nnx 2+nnx;
              2 1 1 1 1 3+2*nnx 1+2*nnx  2+nnx 2+2*nnx 1+nnx];
  Elements = reprow(Elements,1:2,nelx-1,[2 0 0 0 2 2 2 2 2 2]);
  Elements = reprow(Elements,1:2*nelx,nely-1,[2*nelx 0 0 0 2*nnx 2*nnx 2*nnx 2*nnx 2*nnx 2*nnx]);

  % BOUNDARY CONDITIONS AND LOADS
  DOF = getdof(Elements,Types);
  DOF = removedof(DOF,NodesL(:,1));
  F = full(nodalvalues(DOF,NodesB(end,1)+0.02,P));

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
