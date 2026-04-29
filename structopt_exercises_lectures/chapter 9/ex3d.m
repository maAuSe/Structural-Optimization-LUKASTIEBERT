%EX3D
%
%   Topology optimization of an MBB beam for minimum compliance using density
%   filtering in combination with IPOPT (https://github.com/ebertolazzi/mexIPOPT).

% Mattias Schevenels
% January 2024

function ex3d

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

  % RUN OPTIMIZATION
  funcs = struct;
  funcs.objective = @objective;
  funcs.gradient = @gradient;
  funcs.constraints = @constraints;
  funcs.jacobian = @jacobian;
  funcs.jacobianstructure = @jacobianstructure;
  options = struct;
  options.lb = xmin;
  options.ub = xmax;
  options.cl = -inf;
  options.cu = 0;
  options.ipopt.hessian_approximation = 'limited-memory';
  options.ipopt.mu_strategy = 'adaptive';         % See Rojas-Labanda benchmarking paper
  options.ipopt.limited_memory_max_history = 25;  % See Rojas-Labanda benchmarking paper
  options.ipopt.nlp_scaling_method = 'none';      % See Rojas-Labanda benchmarking paper
  options.ipopt.alpha_for_y = 'full';             % See Rojas-Labanda benchmarking paper
  options.ipopt.recalc_y = 'yes';                 % See Rojas-Labanda benchmarking paper
  options.ipopt.tol = 1e-4;
  options.ipopt.print_level = 0;
  xold = nan(nely,nelx);            % Define as global variable
  xcache = nan(nely,nelx);          % Define as global variable
  xPhys = nan(nely,nelx);           % Define as global variable
  ce = nan(nely,nelx);              % Define as global variable
  c = nan;                          % Define as global variable
  iter = -2;                        % IPOPT evaluates gradient twice before iterating
  funevals = 0;
  [x, info] = ipopt(x,funcs,options);

  % COPY VARIABLES TO CALLER WORKSPACE
  S = whos;
  for k = 1:length(S)
    assignin('caller',S(k).name,eval(S(k).name));
  end

%-------------------------------------------------------------------------------

  % OBJECTIVE FUNCTION
  function f0 = objective(x)
    funevals = funevals+1;
    xcache = x;

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
  end

  function df0dx = gradient(x)
    iter = iter+1;

    % SOLVE STATE PROBLEM IF NEEDED
    if ~all(x==xcache)
      objective(x);
    end

    % SENSITIVITIES
    dcdxPhys = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dcdx = conv2(dcdxPhys./Hs,h,'same');
    df0dx = dcdx(:)';

    % PLOT DESIGN
    colormap(gray);
    imagesc(1-xPhys);
    caxis([0 1]);
    axis('equal');
    axis('off');
    drawnow;

    % PRINT STATUS INFORMATION
    change = max(abs(x(:)-xold(:)));
    xold = x;
    if iter>=0
      fprintf('[%s]  iter: %5i | change: % 11.5f |  funevals: %5i | f0: % 11.5f\n',datestr(now,'HH:MM:SS'),iter,change,funevals,c);
    end
  end


%-------------------------------------------------------------------------------

  % CONSTRAINT
  function f = constraints(x)

    % APPLY FILTER TO OBTAIN PHYSICAL DENSITIES
    xPhys = conv2(x,h,'same')./Hs;

    % CONSTRAINT
    v = sum(xPhys(:))/nelx/nely;
    f = v/volfrac-1;
  end

  function dfdx = jacobian(x)

    % CONSTRAINT SENSITIVITIES
    dvdxPhys = repmat(1/nelx/nely,nely,nelx);
    dvdx = conv2(dvdxPhys./Hs,h,'same');
    dfdx = dvdx(:)'/volfrac;
    dfdx = sparse(dfdx);
  end

  function J = jacobianstructure()
    J = sparse(ones(1,nelx*nely));
  end

%-------------------------------------------------------------------------------

end