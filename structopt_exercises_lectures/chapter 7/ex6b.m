%EX6B
%
%   Size and shape optimization of a gridshell for minimum compliance subject to
%   a volume constraint - spline representation, with an extra constraint on the
%   entrance height.

%   This code contains the solution of an exercise introduced in the web lectures.
%   Try to solve the exercise independently before inspecting or running the code.

% Mattias Schevenels
% February 2021

% QUARTER GRIDSHELL DIMENSIONS
Lx = 5;                                                     % Width
Ly = 5;                                                     % Depth
nx = 8;                                                     % Number of subdivisions in x-direction (must be even)
ny = 8;                                                     % Number of subdivisions in y-direction (must be even)
H = 15;                                                     % Maximum vertex height

% DESIGN VARIABLES
ns = 8;                                                     % Number of shape parameters xs (assuming that Nu=Nv=5)
nh = 1;                                                     % Number of size parameters xh
n = ns+nh;                                                  % Total number of design variables
xs = repmat(H/3,ns,1);                                      % Initial values shape parameters
xsmin = zeros(ns,1);                                        % Minimum values shape parameters
xsmax = repmat(H,ns,1);                                     % Maximum values shape parameters
xh = 0.04;                                                  % Initial value size parameter
xhmin = 0.008;                                              % Minimum value size parameter
xhmax = 0.2;                                                % Maximum value size parameter
x = [xs; xh];                                               % Initial values design variables
xmin = [xsmin; xhmin];                                      % Minimum values design variables
xmax = [xsmax; xhmax];                                      % Maximum values design variables

% PARAMETRIC SHAPE REPRESENTATION
p = 4;                                                      % Spline degree
u = linspace(0.5,1,nx+1);                                   % Free parameter values corresponding to FE nodes in first dimension
v = linspace(0.5,1,ny+1);                                   % Free parameter values corresponding to FE nodes in second dimension
Bu = bspline(p,[],u);                                       % Spline basis functions in first dimension
Bv = bspline(p,[],v);                                       % Spline basis functions in second dimension
Nu = size(Bu,2);                                            % Number of control vertices along first dimension
Nv = size(Bv,2);                                            % Number of control vertices along second dimension
N = Nu*Nv;                                                  % Total number of control vertices
AsMap = [0 3 6 3 0                                          % Map spline control vertices to design variables (assuming than Nu=Nv=5)
         1 4 7 4 1
         2 5 8 5 2
         1 4 7 4 1
         0 3 6 3 0];
iAs = 1:numel(AsMap);                                       % Row indices of non-zeros in As (used for zs = As*x)
jAs = max(1,AsMap(:));                                      % Column indices of non-zeros in As
sAs = min(1,AsMap(:));                                      % Values of non-zeros in As
As = sparse(iAs,jAs,sAs,N,n);                               % Coefficient matrix As
zs = As*x;                                                  % Control vertex z-coordinates
Zs = reshape(zs,Nu,Nv);                                     % Idem in matrix format
Zn = Bu*Zs*Bv';                                             % Nodal z-coordinates in matrix format
zn = Zn(:);                                                 % Idem in vector format

% PARAMETRIC SIZE REPRESENTATION
Ah = zeros(1,n);                                            % Coefficient matrix Ah (used for h = Ah*x)
Ah(end) = 1;
h = Ah*x;                                                   % Bar section side length

% CONSTRAINTS
Vmax = 0.2;                                                 % Maximum material volume
zBmin = 5;                                                  % Minimum height node B
zCmin = 5;                                                  % Minimum height node C

% NODES [iNode x y z]
[Xn,Yn] = ndgrid(linspace(0,Lx,nx+1),linspace(0,Ly,ny+1));
xn = Xn(:);
yn = Yn(:);
nNode = numel(xn)+1;
Nodes = [[1:nNode]' [xn; Lx/2] [yn; Ly/2] [zn; Lx*1e3]];

% ELEMENT TYPES {iType type}
Types = {1 'beam';
         2 'truss'};

% SECTIONS [iSec A ky kz Ixx Iyy Izz]
beta = 0.1406;
Sections = [1 h^2   inf inf beta*h^4   h^4/12 h^4/12;
            2 h^2/2 inf inf beta*h^4/2 h^4/24 h^4/24];

% MATERIALS [iMat E nu]
Materials = [1 10e9 50];                                    % For Euler beams, nu is only used to compute G to determine torsion stiffness

% ELEMENTS [iElt iType iSec iMat iNode1 iNode2 iNode3]
Elements = [1 1 1 1 1    2    nNode;
            2 1 1 1 1    nx+2 nNode;
            3 2 1 1 1    nx+3 nNode;
            4 1 1 1 2    3    nNode;
            5 1 1 1 2    nx+3 nNode;
            6 2 1 1 nx+3 3    nNode];
Elements = reprow(Elements,1:6,nx/2-1,[6 0 0 0 2 2 0]);
Elements = [Elements; 3*nx+1 1 1 1 nx+1 2*nx+2 nNode];
Elements = [Elements;
            3*nx+2 1 1 1 nx+2   nx+3   nNode;
            3*nx+3 1 1 1 nx+2   2*nx+3 nNode;
            3*nx+4 2 1 1 2*nx+3 nx+3   nNode;
            3*nx+5 1 1 1 nx+3   nx+4   nNode;
            3*nx+6 1 1 1 nx+3   2*nx+4 nNode;
            3*nx+7 2 1 1 nx+3   2*nx+5 nNode];
Elements = reprow(Elements,3*nx+2:3*nx+7,nx/2-1,[6 0 0 0 2 2 0]);
Elements = [Elements; 6*nx+2 1 1 1 2*nx+2 3*nx+3 nNode];
Elements = reprow(Elements,1:6*nx+2,ny/2-1,[6*nx+2 0 0 0 2*nx+2 2*nx+2 0]);
Elements = [Elements; (3*nx+1)*ny+1 1 1 1 ny*(nx+1)+1 ny*(nx+1)+2 nNode];
Elements = reprow(Elements,(3*nx+1)*ny+1,nx-1,[1 0 0 0 1 1 0]);
xElem = (Nodes(Elements(:,5),2)+Nodes(Elements(:,6),2))/2;
yElem = (Nodes(Elements(:,5),3)+Nodes(Elements(:,6),3))/2;
zElem = (Nodes(Elements(:,5),4)+Nodes(Elements(:,6),4))/2;
Elements(abs(xElem-0)<1e-8,3) = 2;
Elements(abs(yElem-0)<1e-8,3) = 2;

% DEGREES OF FREEDOM
DOF = getdof(Elements,Types);
[~,jNode] = ismember(floor(DOF),Nodes(:,1));
xDOF = Nodes(jNode,2);
yDOF = Nodes(jNode,3);
zDOF = Nodes(jNode,4);
dDOF = round(100*(DOF-floor(DOF)));
iFixedDOF = find(((xDOF==Lx)&(yDOF==Ly)&(dDOF<=3)) | ...
                 ((xDOF==0)&((dDOF==1)|(dDOF==5)|(dDOF==6))) | ...
                 ((yDOF==0)&((dDOF==2)|(dDOF==4)|(dDOF==6))));
iAllDOF = 1:length(DOF);
iFreeDOF = setdiff(iAllDOF,iFixedDOF);
DOF = DOF(iFreeDOF);
xDOF = xDOF(iFreeDOF);
yDOF = yDOF(iFreeDOF);
zDOF = zDOF(iFreeDOF);
dDOF = dDOF(iFreeDOF);
nDOF = length(DOF);

% LOAD (PHYSICAL AND ADJOINT)
F = zeros(nDOF,1);
F(dDOF==3) = -1000*Lx/nx*Ly/ny;
F(xDOF==0) = F(xDOF==0)/2;
F(yDOF==0) = F(yDOF==0)/2;

% OPTIMIZATION LOOP
iter = 1;                         % Iteration counter
itermax = 1000;                   % Maximum number of iterations
change = inf;                     % (Fictitious) initial design change
tol = 1e-4;                       % Convergence threshold
mmaparams = [];                   % Internal MMA parameters
while true

  % UPDATE NODES MATRIX
  zs = As*x;                                                % Control vertex z-coordinates
  Zs = reshape(zs,Nu,Nv);                                   % Idem in matrix format
  Zn = Bu*Zs*Bv';                                           % Nodal z-coordinates in matrix format
  zn = Zn(:);                                               % Idem in vector format
  Nodes(1:nNode-1,4) = zn;                                  % Nodes matrix

  % UPDATE SECTIONS MATRIX
  h = Ah*x;                                                 % Bar section side length
  Sections(:,[2 5 6 7]) = [h^2   beta*h^4   h^4/12 h^4/12;
                           h^2/2 beta*h^4/2 h^4/24 h^4/24];

  % OBJECTIVE FUNCTION AND CONSTRAINTS
  K = asmkm(Nodes,Elements,Types,Sections,Materials,DOF);
  U = K\F;
  C = F'*U;
  Ve = elemvolumes(Nodes,Elements,Types,Sections);
  V = sum(Ve);
  zA = Nodes(1,4);                                          % z-coordinate of node A
  zB = Nodes(nx+1,4);                                       % z-coordinate of node B
  zC = Nodes(ny*(nx+1)+1,4);                                % z-coordinate of node C
  f0 = C;
  f = [V/Vmax-1; 1-zB/zBmin; 1-zC/zCmin];

  % SENSITIVITIES
  if iter==1 || (f0<=subp.f0e && all(f<=subp.fe))
    dzndzs = kron(Bv,Bu);                                     % Derivative of zn=Zn(:) wrt zs=Zs(:)
    dzsdx = As;                                               % Derivative of zs wrt x
    dzndx = dzndzs*dzsdx;                                     % Derivative of zn wrt x (via chain rule)
    dNodesdx = zeros([size(Nodes),n]);
    dNodesdx(1:nNode-1,4,:) = dzndx;
    dSectionsdx = zeros([size(Sections),n]);
    dSectionsdx(1,2,:) = 2*h*Ah;
    dSectionsdx(1,5,:) = 4*beta*h^3*Ah;
    dSectionsdx(1,6,:) = h^3/3*Ah;
    dSectionsdx(1,7,:) = h^3/3*Ah;
    dSectionsdx(2,2,:) = 2*h*Ah/2;
    dSectionsdx(2,5,:) = 4*beta*h^3*Ah/2;
    dSectionsdx(2,6,:) = h^3/3*Ah/2;
    dSectionsdx(2,7,:) = h^3/3*Ah/2;
    [~,~,dKdx] = asmkm(Nodes,Elements,Types,Sections,Materials,DOF,dNodesdx,dSectionsdx);
    dCdx = zeros(1,n);
    for k = 1:n
      dCdx(k) = -U'*dKdx{k}*U;
    end
    [~,dVedx] = elemvolumes(Nodes,Elements,Types,Sections,dNodesdx,dSectionsdx);
    dVdx = sum(dVedx);
    dzAdx = reshape(dNodesdx(1,4,:),1,[]);
    dzBdx = reshape(dNodesdx(nx+1,4,:),1,[]);
    dzCdx = reshape(dNodesdx(ny*(nx+1)+1,4,:),1,[]);
    df0dx = dCdx;
    dfdx = [dVdx/Vmax; -dzBdx/zBmin; -dzCdx/zCmin];
  end

  % GCMMA UPDATE
  [xnew,y,z,lmult,mmaparams,subp,change,history] = gcmma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams);

  % PLOT CURRENT DESIGN
  if ~isinf(change)
    Nodes2 = Nodes; Nodes2(:,2) = -Nodes2(:,2);
    Nodes3 = Nodes; Nodes3(:,3) = -Nodes3(:,3);
    Nodes4 = Nodes; Nodes4(:,2:3) = -Nodes4(:,2:3);
    plotelem(Nodes,Elements,Types,'Numbering','off','GCS','off');
    hold('on');
    plotelem(Nodes2,Elements,Types,'Numbering','off','GCS','off');
    plotelem(Nodes3,Elements,Types,'Numbering','off','GCS','off');
    plotelem(Nodes4,Elements,Types,'Numbering','off','GCS','off');
    hold('off');
    set(gca,'Position',[0 0 1 1]);
    axis([-5 5 -5 5 -1 7]);
    drawnow;
  end

  % CHECK CONVERGENCE CRITERIA AND MOVE ON TO NEXT ITERATION
  if change>tol && (iter<itermax || isinf(change))
    iter = iter+1;
    x = xnew;
  else
    break;
  end
end

