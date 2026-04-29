function [xnew,xgbest,f0gbest,fgbest,psoparams,change,history] = pso(x,xmin,xmax,f0,f,psoparams,silent)

%PSOSUB   Particle swarm optimization.
%
%   [xnew,xgbest,f0gbest,fgbest,psoparams,change,history] = pso(x,xmin,xmax,f0,f,psoparams)
%   performs one particle swarm optimization iteration.  Constraints are taken
%   into account by means of extremely strong penalization.
%
%   INPUT ARGUMENTS
%
%   x                 Current particle locations (n * p).
%   xmin              Lower bounds for x (n * 1) or (1 * 1).
%   xmax              Upper bounds for x (n * 1) or (1 * 1).
%   f0                Value of the objective function at x (1 * p).
%   f                 Values of the constraints at x (m * p).
%   psoparams         Internal parameters controlling the evolution of the particle
%                     swarm.  Use [] in the first iteration and return the
%                     corresponding output argument in the following iterations.
%
%   OUTPUT ARGUMENTS
%
%   xnew              Updated particle locations (n * p).
%   xgbest            Updated global best location (n * 1).
%   f0gbest           Value of the objective function at xgbest (1 * 1).
%   fgbest            Values of the constraints at xgbest (m * 1).
%   psoparams         Updated values of the internal PSO parameters (1 * 3).
%   change            Design change, defined as MAX(MAX(ABS(xnew-x))).
%   history.time      Time at iteration i (1 * N).
%   history.iter      Iteration counter (1 * N).
%   history.funevals  Total number of function evaluations performed at iteration i (1 * N).
%   history.x         Values of variables x at iteration i (n * p * N).
%                     If the n * p exceeds 1000, an empty matrix is returned.
%   history.change    Design change, defined as above (1 * N).
%   history.xgbest    Global best location at iteration i (n * N).
%                     If the number of design variables exceeds 100, an empty matrix is returned.
%   history.f0gbest   Value of objective function at xgbest in iteration i (1 * N).
%   history.fgbest    Value of constraints at xgbest in iteration i (m * N).
%                     If the number of constraints exceeds 100, MAX(f) is returned instead (1 * N).
%
%   The following information is printed as output to the screen:
%
%   iter:             Iteration counter.
%   funevals:         Total number of function evaluations performed.
%   change:           Design change, defined as above.
%   f0:               Value of the objective function at xgbest.
%   fmax:             Value of the most critical constraint at xgbest, i.e. MAX(fgbest).
%
%   Use PSO(...,'silent') to avoid printing output.

% Mattias Schevenels
% February 2021

% This function is partly based on the file particleswarm.m from the MATLAB
% Global Optimization Toolbox.

% CHECK SILENT STATUS
silent = exist('silent','var') && strcmpi(silent,'silent');

% NUMBER OF CONSTRAINTS, DESIGN VARIABLES, AND PARTICLES
m = size(f,1);                 % Number of constraints
n = size(x,1);                 % Number of design variables
p = size(x,2);                 % Population size

% VECTORIZE SCALAR INPUT ARGUMENTS
if numel(xmin)==1, xmin = repmat(xmin,n,1); end
if numel(xmax)==1, xmax = repmat(xmax,n,1); end

% TRANSPOSE SINGLE-ROW INPUT ARGUMENTS
if size(xmin,1)==1, xmin = xmin.'; end
if size(xmax,1)==1, xmax = xmax.'; end
if size(f0,2)==1; f0 = f0.'; end

% CHECK INPUT ARGUMENT DIMENSIONS
checkdim(x,[n p],'Input argument x must have dimensions (n * p).');
checkdim(xmin,[n 1],'Input argument xmin must have dimensions (1 * 1) or (n * 1) where n = SIZE(x,1).');
checkdim(xmax,[n 1],'Input argument xmax must have dimensions (1 * 1) or (n * 1) where n = SIZE(x,1).');
checkdim(f0,[1 p],'Input argument f0 must have dimensions (1 * p) where p = SIZE(x,2).');
checkdim(f,[m p],'Input argument f must have dimensions (m * p) where p = SIZE(x,2).');

% FIXED PSO PARAMETERS
y1 = 1.49;                     % Self adjustment weight
y2 = 1.49;                     % Social adjustment weight
Wmin = 0.1;                    % Minimum inertia
Wmax = 1.1;                    % Maximum inertia
Nmin = max(1,floor(0.25*p));   % Minimum neighborhood size

% INITIALIZE INTERNAL PSO PARAMETERS ON FIRST ITERATION
if isempty(psoparams)
  psoparams = struct;
  psoparams.iter = 0;
  psoparams.C = 0;
  psoparams.W = Wmax;
  psoparams.N = Nmin;
  psoparams.v = repmat(xmin-xmax,1,p)+2*rand(n,p).*repmat(xmax-xmin,1,p);
  psoparams.xpbest = nan(n,p);
  psoparams.f0pbest = inf(1,p);
  psoparams.fpbest = inf(m,p);
  psoparams.history.time = [];
  psoparams.history.iter = [];
  psoparams.history.funevals = [];
  psoparams.history.x = [];
  psoparams.history.change = [];
  psoparams.history.xgbest = [];
  psoparams.history.f0gbest = [];
  psoparams.history.fgbest = [];
end

% RETRIEVE PSO PARAMETERS FROM PREVIOUS ITERATION
iter = psoparams.iter;
C = psoparams.C;
W = psoparams.W;
N = psoparams.N;
v = psoparams.v;
xpbest = psoparams.xpbest;
f0pbest = psoparams.f0pbest;
fpbest = psoparams.fpbest;
history = psoparams.history;

% INCREASE ITERATION COUNTER
iter = iter+1;

% DETERMINE PENALIZED FITNESS
alpha = sqrt(realmax);
F = f0+alpha*sum(max(0,f).^2,1);
Fpbest = f0pbest+alpha*sum(max(0,fpbest).^2,1);

% UPDATE EACH PARTICLE'S PERSONAL BEST
flag = min(F)<min(Fpbest);
k = F<Fpbest;
xpbest(:,k) = x(:,k);
f0pbest(k) = f0(k);
if m>0, fpbest(:,k) = f(:,k); end
Fpbest(k) = F(k);

% UPDATE GLOBAL BEST
[~,k] = min(Fpbest);
xgbest = xpbest(:,k);
f0gbest = f0pbest(k);
fgbest = fpbest(:,k);
Fgbest = Fpbest(k);

% UPDATE NEIGHBORHOOD PARAMETERS
if flag && iter>1
  C = max(0,C-1);
  N = Nmin;
elseif iter>1
  C = C+1;
  N = min(p,N+Nmin);
end
if C<2, W = 2*W; end
if C>5, W = W/2; end
W = max(Wmin,min(Wmax,W));

% DETERMINE EACH PARTICLE'S NEIGHBORHOOD BEST
neighborIndex = zeros(N,p);
neighborIndex(1,:) = 1:p;
for k = 1:p
  neighbors = randperm(p-1, N-1);
  iShift = neighbors >= k;
  neighbors(iShift) = neighbors(iShift) + 1;
  neighborIndex(2:end,k) = neighbors;
end
[~,k] = min(Fpbest(neighborIndex));
k = k+[0:N:(p-1)*N];
k = neighborIndex(k);
xnbest = xpbest(:,k);

% UPDATE PARTICLE VELOCITIES AND POSITIONS
% v(:,F>f0) = 0;
v = W*v + y1*rand(n,p).*(xpbest-x) + y2*rand(n,p).*(xnbest-x);
xnew = x+v;
xnew = max(xnew,repmat(xmin,1,p));
xnew = min(xnew,repmat(xmax,1,p));
v(xnew==repmat(xmin,1,p)) = 0;
v(xnew==repmat(xmax,1,p)) = 0;
change = max(max(abs(xnew-x)));

% PRINT CONVERGENCE HISTORY INFORMATION
time = now;
if ~silent
  fprintf('[%s]  iter: %5i | funevals: %7i | change: % 10.5f | f0: % 10.5f | fmax: % 10.5f\n',datestr(time,'HH:MM:SS'),iter,p*iter,change,f0gbest,max(fgbest));
end

% COLLECT CONVERGENCE HISTORY INFORMATION
if nargout>=7
  history.time(:,iter) = time;
  history.iter(:,iter) = iter;
  history.funevals(:,iter) = p*iter;
  if n*p<=1000, history.x(:,:,iter) = x; end
  history.change(:,iter) = change;
  if n<=100, history.xgbest(:,iter) = xgbest; end
  history.f0gbest(:,iter) = f0gbest;
  if m<=100, history.fgbest(:,iter) = fgbest; else history.fgbest(:,iter) = max(fgbest); end
end

% COLLECT INTERNAL PSO PARAMETERS
psoparams = struct;
psoparams.iter = iter;
psoparams.C = C;
psoparams.W = W;
psoparams.N = N;
psoparams.v = v;
psoparams.xpbest = xpbest;
psoparams.f0pbest = f0pbest;
psoparams.fpbest = fpbest;
psoparams.history = history;

% SUBFUNCTION CHECKDIM - CHECK INPUT ARGUMENT DIMENSIONS
function checkdim(x,dim,err)
if length(size(x))~=length(dim) || any(size(x)~=dim), error(err); end


