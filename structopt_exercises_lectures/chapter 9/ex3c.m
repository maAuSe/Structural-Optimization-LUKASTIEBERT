%EX3C
%
%   Topology optimization of an MBB beam for minimum compliance using density
%   filtering in combination with FMINCON (interior point method).

% Mattias Schevenels
% January 2024

function ex3c

  % PROBLEM SETUP
  nelx = 120;          % Number of elements in x-direction
  nely = 40;           % Number of elements in y-direction
  volfrac = 0.5;       % Maximum volume fraction
  penal = 3;           % Penalization power
  rmin = 4.8;          % Filter radius in terms of elements

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
  xmin = repmat(0,nely,nelx);       % Minimum value design variables
  xmax = repmat(1,nely,nelx);       % Maximum value design variables
  xold = nan(nely,nelx);

  % RUN OPTIMIZATION
  options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'StepTolerance', 0.01, ...
    'SubproblemAlgorithm','cg', ...
    'MaxProjCGIter',100, ...
    'HessianApproximation','lbfgs', ...
    'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true, ...
    'OutputFcn',@outfun, ...
    'Display','none');
  [x(:),fval,exitflag,output,lambda,grad] = fmincon(@objfun,x,[],[],[],[],xmin,xmax,@nonlcon,options);

  % COPY VARIABLES TO CALLER WORKSPACE
  S = whos;
  for k = 1:length(S)
    assignin('caller',S(k).name,eval(S(k).name));
  end

%-------------------------------------------------------------------------------

  % OBJECTIVE FUNCTION
  function [f0,df0dx] = objfun(x)

    % APPLY FILTER TO OBTAIN PHYSICAL DENSITIES
    xPhys = conv2(x,h,'same')./Hs;

    % FINITE ELEMENT ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);

    % OBJECTIVE FUNCTION
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
    f0 = c;

    % SENSITIVITIES
    dcdxPhys = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;        % Sensitivity wrt physical densities
    dcdx = conv2(dcdxPhys./Hs,h,'same');                     % Sensitivity wrt design variables
    df0dx = dcdx(:);
  end

%-------------------------------------------------------------------------------

  % CONSTRAINT
  function [f,dum1,dfdx,dum2] = nonlcon(x)

    % APPLY FILTER TO OBTAIN PHYSICAL DENSITIES
    xPhys = conv2(x,h,'same')./Hs;

    % CONSTRAINT
    v = sum(xPhys(:))/nelx/nely;
    f = v/volfrac-1;

    % SENSITIVITIES
    dvdxPhys = repmat(1/nelx/nely,nely,nelx);                % Sensitivity wrt physical densities
    dvdx = conv2(dvdxPhys./Hs,h,'same');                     % Sensitivity wrt design variables
    dfdx = dvdx(:)/volfrac;

    dum1 = [];
    dum2 = [];
  end

%-------------------------------------------------------------------------------

  % PRINT AND PLOT FUNCTION
  function stop = outfun(x,optimValues,state)

    if strcmpi(state,'iter')

      % PLOT DESIGN
      xPhys = conv2(x,h,'same')./Hs;
      colormap(gray);
      imagesc(1-xPhys);
      caxis([0 1]);
      axis('equal');
      axis('off');
      drawnow;

      % PRINT STATUS INFORMATION
      change = max(abs(x(:)-xold(:)));
      xold = x;
      fprintf('[%s]  iter: %5i | change: % 11.5f |  funevals: %5i | f0: % 11.5f | constrviolation: % 11.5f\n',datestr(now,'HH:MM:SS'),optimValues.iteration,change,optimValues.funccount,optimValues.fval,optimValues.constrviolation);
    end

    % DON'T STOP
    stop = false;
  end

%-------------------------------------------------------------------------------

end