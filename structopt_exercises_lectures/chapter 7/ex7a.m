%EX7A
%
%   Shape optimization of a plane sheet for minimum compliance, using a
%   B-spline surface mesh.

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
pu = 2;          % Spline degree in terms of u (horizontal direction)
pv = 1;          % Spline degree in terms of v (vertical direction)
Nu = 6;          % Number of control points in u-direction
Nv = 2;          % Number of control points in v-direction
n = Nu;          % Number of design variables (y-coordinates of top control points)

% DESIGN VARIABLES (y-coordinates of control points along top boundary, ordered left -> right)
x = repmat(H,n,1); 
xmin = H/5;                               
xmax = H*2;

% SPLINE SAMPLING POINTS (top+bottom / left+right boundary)
u = linspace(0,1,nnx);
v = linspace(0,1,nny);

% SPLINE BASIS FUNCTIONS 
Bu = bspline(pu,Nu,u);          
Bv = bspline(pv,Nv,v);          

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

  % CONTROL POINT COORDINATES (top / left / bottom / right boundary)
  [X,Y] = ndgrid(linspace(0,W,Nu),linspace(0,H,Nv));
  Y(:,end) = x;
  
  % CONTROL POINT COORDINATE SENSITIVITIES
  dXdx = zeros(Nu*Nv,n);
  dYdx = zeros(Nu*Nv,n);
  dYdx(Nu*(Nv-1)+1:end,1:n) = eye(n);

  % NODE COORDINATES
  xNode = Bu*X*Bv';  xNode = xNode(:);
  yNode = Bu*Y*Bv';  yNode = yNode(:);
  Nodes = [[1:nnx*nny]' xNode yNode zeros(nnx*nny,1)];
  
  % NODE COORDINATE SENSITIVITIES
  dxNodedX = kron(Bv,Bu);
  dyNodedY = kron(Bv,Bu);
  dxNodedx = dxNodedX*dXdx;
  dyNodedx = dyNodedY*dYdx;
  dNodesdx = zeros(nnx*nny,4,n);
  dNodesdx(:,2,:) = reshape(dxNodedx,[],1,n);
  dNodesdx(:,3,:) = reshape(dyNodedx,[],1,n);

  % ELEMENTS
  Elements = [1 1 1 1 1 3 3+2*nnx  2 3+nnx 2+nnx;
              2 1 1 1 1 3+2*nnx 1+2*nnx  2+nnx 2+2*nnx 1+nnx];
  Elements = reprow(Elements,1:2,nelx-1,[2 0 0 0 2 2 2 2 2 2]);
  Elements = reprow(Elements,1:2*nelx,nely-1,[2*nelx 0 0 0 2*nnx 2*nnx 2*nnx 2*nnx 2*nnx 2*nnx]);

  % BOUNDARY CONDITIONS AND LOADS
  DOF = getdof(Elements,Types);
  DOF = removedof(DOF,1:nnx:nnx*nny);
  F = full(nodalvalues(DOF,nnx+0.02,P));

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
  plot(X,Y,'ro','MarkerFaceColor','r','MarkerSize',4);
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
