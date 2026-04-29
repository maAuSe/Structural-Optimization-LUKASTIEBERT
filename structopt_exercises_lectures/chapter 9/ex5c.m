%EX5C
%
%   Topology optimization considering different load cases - the maximum of the
%   compliance values is minimized.

%   This code contains the solution of an exercise introduced in the web lectures.
%   Try to solve the exercise independently before inspecting or running the code.

% Mattias Schevenels
% February 2021

% PROBLEM SETUP
nelx = 100;          % Number of elements in x-direction
nely = 100;          % Number of elements in y-direction
volfrac = 0.2;       % Maximum volume fraction
penal = 3;           % Penalization power
rmin = 3;            % Filter radius in terms of elements

% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;

% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

% DEFINE LOADS AND SUPPORTS
F1 = sparse(2,1,1,2*(nely+1)*(nelx+1),1);
F2 = sparse(2*(nely+1),1,-1,2*(nely+1)*(nelx+1),1);
F = [F1 F2];
U = zeros(2*(nely+1)*(nelx+1),2);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

% PREPARE FILTER
[dy,dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
h = max(0,rmin-sqrt(dx.^2+dy.^2));
Hs = conv2(ones(nely,nelx),h,'same');

% DESIGN VARIABLES
x = repmat(volfrac,nely,nelx);    % Initialization of design variables
xmin = 0;                         % Minimum value design variables
xmax = 1;                         % Maximum value design variables

% MMA PARAMETERS
a0 = 1;
a = [1;1;0];
c = 10000;
d = 0;

% OPTIMIZATION LOOP
iter = 1;                         % Iteration counter
change = 1;                       % (Fictitious) initial design change
tol = 0.01;                       % Convergence threshold
mmaparams = [];                   % Internal MMA parameters
while change>tol

  % APPLY FILTER TO OBTAIN PHYSICAL DENSITIES
  xPhys = conv2(x,h,'same')./Hs;

  % FINITE ELEMENT ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
  U1 = U(:,1);
  U2 = U(:,2);

  % OBJECTIVE FUNCTION AND CONSTRAINT
  ce1 = reshape(sum((U1(edofMat)*KE).*U1(edofMat),2),nely,nelx);
  ce2 = reshape(sum((U2(edofMat)*KE).*U2(edofMat),2),nely,nelx);
  c1 = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce1));
  c2 = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce2));
  v = sum(xPhys(:))/nelx/nely;
  f0 = 0;
  f = zeros(3,1);
  f(1) = c1/100;
  f(2) = c2/100;
  f(3) = v/volfrac-1;

  % SENSITIVITIES
  dc1dxPhys = -penal*(E0-Emin)*xPhys.^(penal-1).*ce1;      % Sensitivity wrt physical densities
  dc2dxPhys = -penal*(E0-Emin)*xPhys.^(penal-1).*ce2;      % Sensitivity wrt physical densities
  dvdxPhys = repmat(1/nelx/nely,nely,nelx);                % Sensitivity wrt physical densities
  dc1dx = conv2(dc1dxPhys./Hs,h,'same');                   % Sensitivity wrt design variables
  dc2dx = conv2(dc2dxPhys./Hs,h,'same');                   % Sensitivity wrt design variables
  dvdx = conv2(dvdxPhys./Hs,h,'same');                     % Sensitivity wrt design variables
  df0dx = zeros(nelx*nely,1);
  dfdx = zeros(3,nelx*nely);
  dfdx(1,:) = dc1dx(:)'/100;
  dfdx(2,:) = dc2dx(:)'/100;
  dfdx(3,:) = dvdx(:)'/volfrac;
  % MMA UPDATE
  [xnew,~,~,~,mmaparams,~,change,history] = mma(x(:),xmin,xmax,f0,f,df0dx,dfdx,mmaparams,[],[],[],[],[],[],[],a0,a,c,d,'silent');

  % PRINT CONVERGENCE HISTORY INFORMATION
  fprintf('[%s]  iter: %5i | change: % 11.5f | c1: % 11.5f | c2: % 11.5f | v: % 11.5f\n',datestr(now,'HH:MM:SS'),iter,change,c1,c2,v);

  % PLOT CURRENT DESIGN
  colormap(gray);
  imagesc(1-xPhys);
  caxis([0 1]);
  axis('equal');
  axis('off');
  drawnow;

  % CHECK CONVERGENCE CRITERIA AND MOVE ON TO NEXT ITERATION
  if change>tol
    iter = iter+1;
    x(:) = xnew;
  end

end
