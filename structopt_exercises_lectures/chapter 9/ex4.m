%EX4A
%
%   Topology optimization of an MBB beam for minimum compliance using
%   sensitivity filtering, density filtering, or Heaviside projection filtering,
%   in combination with MMA.

% Mattias Schevenels
% February 2021

% PROBLEM SETUP
nelx = 120;          % Number of elements in x-direction
nely = 40;           % Number of elements in y-direction
volfrac = 0.5;       % Maximum volume fraction
penal = 3;           % Penalization power
rmin = 4.8;          % Filter radius in terms of elements
ft = 3;              % Filter type: 1 for sensitivity, 2 for density, 3 for Heaviside
eta = 0.5;           % Heaviside projection threshold

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

% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
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

% OPTIMIZATION LOOP
beta = 1;                         % Initial Heaviside projection smoothness parameter value
betamax = 16;                     % Final beta value
iter = 1;                         % Iteration counter
iterbeta = 1;                     % Iteration counter since last increase of beta
iterbetamax = 50;                 % Number of iterations per beta value
change = 1;                       % (Fictitious) initial design change
tol = 0.01;                       % Convergence threshold
mmaparams = [];                   % Internal MMA parameters
while change>tol || (ft == 3 && beta<betamax)

  % APPLY FILTER TO OBTAIN PHYSICAL DENSITIES
  if ft == 1
    xPhys = x;
  elseif ft == 2
    xPhys = conv2(x,h,'same')./Hs;
  elseif ft == 3
    xTilde = conv2(x,h,'same')./Hs;
    xPhys = (tanh(beta*eta)+tanh(beta*(xTilde-eta)))./(tanh(beta*eta)+tanh(beta*(1-eta)));
  end

  % FINITE ELEMENT ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);

  % OBJECTIVE FUNCTION AND CONSTRAINT
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  v = sum(xPhys(:))/nelx/nely;
  f0 = c/100;
  f = v/volfrac-1;

  % SENSITIVITIES
  if ft == 1
    dcdx = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dvdx = repmat(1/nelx/nely,nely,nelx);
    dcdx = conv2(dcdx.*xPhys,h,'same')./Hs./max(1e-3,xPhys);                                 % Apply sensitivity filter
  elseif ft == 2
    dcdxPhys = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;                                        % Sensitivity wrt physical densities
    dvdxPhys = repmat(1/nelx/nely,nely,nelx);                                                % Sensitivity wrt physical densities
    dcdx = conv2(dcdxPhys./Hs,h,'same');                                                     % Sensitivity wrt design variables
    dvdx = conv2(dvdxPhys./Hs,h,'same');                                                     % Sensitivity wrt design variables
  elseif ft == 3
    dcdxPhys = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;                                        % Sensitivity wrt physical densities xPhys
    dvdxPhys = repmat(1/nelx/nely,nely,nelx);                                                % Sensitivity wrt physical densities xPhys
    dxPhysdxTilde = (beta*sech(beta*(xTilde-eta)).^2./(tanh(beta*eta)+tanh(beta*(1-eta))));  % Derivative dxPhys/dxTilde
    dcdxTilde = dcdxPhys.*dxPhysdxTilde;                                                     % Sensitivity wrt filtered variables xTilde
    dvdxTilde = dvdxPhys.*dxPhysdxTilde;                                                     % Sensitivity wrt filtered variables xTilde
    dcdx = conv2(dcdxTilde./Hs,h,'same');                                                    % Sensitivity wrt design variables
    dvdx = conv2(dvdxTilde./Hs,h,'same');                                                    % Sensitivity wrt design variables
  end
  df0dx = dcdx(:)'/100;
  dfdx = dvdx(:)'/volfrac;

  % MMA UPDATE
  [xnew,~,~,~,mmaparams,~,change,history] = mma(x(:),xmin,xmax,f0,f,df0dx,dfdx,mmaparams);

  % PLOT CURRENT DESIGN
  colormap(gray);
  imagesc(1-xPhys);
  caxis([0 1]);
  axis('equal');
  axis('off');
  drawnow;

  % UPDATE HEAVISIDE REGULARIZATION PARAMETER
  if ft == 3 && beta < betamax && iterbeta >= iterbetamax
    beta = 2*beta;
    iterbeta = 0;
    change = 1;
    fprintf('Parameter beta increased to %g.\n',beta);
  end

  % CHECK CONVERGENCE CRITERIA AND MOVE ON TO NEXT ITERATION
  if change>tol || beta<betamax
    iter = iter+1;
    iterbeta = iterbeta+1;
    x(:) = xnew;
  end

end
