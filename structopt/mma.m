function [xnew,y,z,lmult,mmaparams,subp,change,history] = mma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams,move,low,upp,sinit,sincr,sdecr,mu,a0,a,c,d,arg1)

%MMA   Method of moving asymptotes.
%
%   [xnew,y,z,lmult,mmaparams,subp,change,history] = ...
%   MMA(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams,move,low,upp,sinit,sincr,sdecr,mu,a0,a,c,d)
%   performs one MMA-iteration, aimed at solving the nonlinear programming problem:
%
%      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
%    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
%                xmin_j <= x_j <= xmax_j,    j = 1,...,n
%                z >= 0,   y_i >= 0,         i = 1,...,m
%
%   INPUT ARGUMENTS
%
%   x                 Current values of the variables x (n * 1).
%   xmin              Lower bounds for x (n * 1) or (1 * 1).
%   xmax              Upper bounds for x (n * 1) or (1 * 1).
%   f0                Value of the objective function at x (1 * 1).
%   f                 Values of the constraints at x (m * 1).
%   df0dx             Derivatives of the objective function at x (1 * n).
%   dfdx              Derivatives of the constraints at x (m * n).
%   mmaparams         Internal parameters controlling the moving asymptotes.  Use []
%                     in the first iteration and return the corresponding output
%                     argument in the following iterations.
%   move              Move limit (n * 1) or (1 * 1). Default: 0.2.
%   low               User-specified values for the lower asymptotes (n * 1) or (1 * 1). Default: [].
%   upp               User-specified values for the upper asymptotes (n * 1) or (1 * 1). Default: [].
%                     Specify [] or nan to compute asymptotes automatically.
%                     Use low = 0 and upp = INF to use CONLIN.
%                     Use low = -INF and upp = INF to use SLP.
%                     The next 4 input arguments are not used in the case of CONLIN or SLP.
%   sinit             Asymptotes relative initial position (n * 1) or (1 * 1).  Default: 0.5.
%   sincr             Factor used to increase the distance between the asymptotes (n * 1) or (1 * 1). Default: 1.2.
%   sdecr             Factor used to decrease the distance between the asymptotes (n * 1) or (1 * 1). Default: 0.7.
%   mu                Relative distance to keep away from the asymptotes (n * 1) or (1 * 1). Default: 0.1.
%   a0                The constant a_0 in the term a_0*z (1 * 1). Default: 1.
%   a                 The constants a_i in the terms a_i*z (m * 1) or (1 * 1). Default: 0.
%   c                 The constants c_i in the terms c_i*y_i (m * 1) or (1 * 1). Default: 10000.
%   d                 The constants d_i in the terms 0.5*d_i*(y_i)^2 (m * 1) or (1 * 1). Default: 1.
%
%   OUTPUT ARGUMENTS
%
%   xnew              Optimal values of the variables x in the current MMA subproblem (n * 1).
%   y                 Optimal values of the variables y in the current MMA subproblem (n * 1).
%   z                 Optimal values of the variables z in the current MMA subproblem (n * 1).
%   lmult.lambda      Lagrange multipliers for the general MMA constraints (m * 1).
%   lmult.xi          Lagrange multipliers for the constraints alfa_j - x_j <= 0 (n * 1).
%   lmult.eta         Lagrange multipliers for the constraints x_j - beta_j <= 0 (n * 1).
%   lmult.mu          Lagrange multipliers for the constraints -y_i <= 0 (m * 1).
%   lmult.zeta        Lagrange multiplier for the constraint -z <= 0 (1 * 1).
%   lmult.s           Slack variables for the general MMA constraints (m * 1).
%   mmaparams         Updated values of the internal MMA parameters.
%   subp.p0           Coefficients of the upper asymptote terms in the objective function approximation (n * 1).
%   subp.q0           Coefficients of the lower asymptote terms in the objective function approximation (n * 1).
%   subp.r0           Constant term in the objective function approximation (1 * 1).
%   subp.p            Coefficients of the upper asymptote terms in the constraint approximations (m * n).
%   subp.q            Coefficients of the lower asymptote terms in the constraint approximations (m * n).
%   subp.r            Constant term in the constraint approximations (m * 1).
%   subp.low          Values of the lower asymptotes (n * 1).
%   subp.upp          Values of the upper asymptotes (n * 1).
%   subp.alpha        Lower bounds for x used in the current MMA subproblem (n * 1).
%   subp.beta         Upper bounds for x used in the current MMA subproblem (n * 1).
%   change            Design change, defined as MAX(ABS(x-xold)), where x is the design passed as input in
%                     the current iteration, and xold is the design passed as input in the previous iteration.
%   history.time      Time at iteration i (1 * N).
%   history.iter      Iteration counter (1 * N).
%   history.funevals  Total number of function evaluations performed at iteration i (1 * N).
%   history.x         Values of variables x at iteration i (n * N).
%                     If the number of design variables exceeds 100, an empty matrix is returned.
%   history.change    Design change, defined as above (1 * N).
%   history.f0        Value of objective function at iteration i (1 * N).
%   history.f         Value of constraints at iteration i (m * N).
%                     If the number of constraints exceeds 100, MAX(f) is returned instead (1 * N).
%
%   The following information is printed as output to the screen:
%
%   iter:             Iteration counter.
%   change:           Design change, defined as above.
%   f0:               Value of the objective function at x.
%   fmax:             Value of the most critical constraint at x, i.e. MAX(f).
%
%   Use MMA(...,'silent') to avoid printing output.

% Mattias Schevenels
% January 2026

% This function is a wrapper function around Svanberg's implementation, which is
% included below with minor modifications (denoted by the tag MS).

% DETERMINE SILENT STATUS
silent = false;
if exist('move','var') && strcmpi(move,'silent'), silent=true; clear('move'); end
if exist('low','var') && strcmpi(low,'silent'), silent=true; clear('low'); end
if exist('upp','var') && strcmpi(upp,'silent'), silent=true; clear('upp'); end
if exist('sinit','var') && strcmpi(sinit,'silent'), silent=true; clear('sinit'); end
if exist('sincr','var') && strcmpi(sincr,'silent'), silent=true; clear('sincr'); end
if exist('sdecr','var') && strcmpi(sdecr,'silent'), silent=true; clear('sdecr'); end
if exist('mu','var') && strcmpi(mu,'silent'), silent=true; clear('mu'); end
if exist('a0','var') && strcmpi(a0,'silent'), silent=true; clear('a0'); end
if exist('a','var') && strcmpi(a,'silent'), silent=true; clear('a'); end
if exist('c','var') && strcmpi(c,'silent'), silent=true; clear('c'); end
if exist('d','var') && strcmpi(d,'silent'), silent=true; clear('d'); end
if exist('arg1','var') && strcmpi(arg1,'silent'), silent=true; clear('arg1'); end

% DEFAULT INPUT ARGUMENTS
if ~exist('move','var') || isempty(move), move = 0.2; end
if ~exist('low','var') || isempty(low), low = nan; end
if ~exist('upp','var') || isempty(upp), upp = nan; end
if ~exist('sinit','var') || isempty(sinit), sinit = 0.5; end
if ~exist('sincr','var') || isempty(sincr), sincr = 1.2; end
if ~exist('sdecr','var') || isempty(sdecr), sdecr = 0.7; end
if ~exist('mu','var') || isempty(mu), mu = 0.1; end
if ~exist('a0','var') || isempty(a0), a0 = 1; end
if ~exist('a','var') || isempty(a), a = 0; end
if ~exist('c','var') || isempty(c), c = 10000; end
if ~exist('d','var') || isempty(d), d = 1; end

% NUMBER OF CONSTRAINTS AND DESIGN VARIABLES
m = numel(f);
n = numel(x);

% VECTORIZE SCALAR INPUT ARGUMENTS
if numel(xmin)==1, xmin = repmat(xmin,n,1); end
if numel(xmax)==1, xmax = repmat(xmax,n,1); end
if numel(a)==1, a = repmat(a,m,1); end
if numel(c)==1, c = repmat(c,m,1); end
if numel(d)==1, d = repmat(d,m,1); end
if numel(move)==1, move = repmat(move,n,1); end
if numel(low)==1, low = repmat(low,n,1); end
if numel(upp)==1, upp = repmat(upp,n,1); end
if numel(sinit)==1, sinit = repmat(sinit,n,1); end
if numel(sincr)==1, sincr = repmat(sincr,n,1); end
if numel(sdecr)==1, sdecr = repmat(sdecr,n,1); end
if numel(mu)==1, mu = repmat(mu,n,1); end

% TRANSPOSE SINGLE-ROW INPUT ARGUMENTS
if size(x,1)==1, x = x.'; end;
if size(xmin,1)==1, xmin = xmin.'; end
if size(xmax,1)==1, xmax = xmax.'; end
if size(f,1)==1; f = f.'; end
if size(df0dx,1)==1, df0dx = df0dx.'; end
if n==1 && size(dfdx,1)==1, dfdx = dfdx.'; end
if m==1 && size(dfdx,2)==1, dfdx = dfdx.'; end
if size(a,1)==1, a = a.'; end
if size(c,1)==1, c = c.'; end
if size(d,1)==1, d = d.'; end
if size(move,1)==1, move = move.'; end
if size(low,1)==1, low = low.'; end
if size(upp,1)==1, upp = upp.'; end
if size(sinit,1)==1, sinit = sinit.'; end
if size(sincr,1)==1, sincr = sincr.'; end
if size(sdecr,1)==1, sdecr = sdecr.'; end
if size(mu,1)==1, mu = mu.'; end

% CHECK INPUT ARGUMENT DIMENSIONS
checkdim(x,[n 1],'Input argument x must have dimensions (n * 1) where n = NUMEL(x).');
checkdim(xmin,[n 1],'Input argument xmin must have dimensions (1 * 1) or (n * 1) where n = NUMEL(x).');
checkdim(xmax,[n 1],'Input argument xmax must have dimensions (1 * 1) or (n * 1) where n = NUMEL(x).');
checkdim(f0,[1 1],'Input argument f0 must have dimensions (1 * 1).');
if m>0, checkdim(f,[m 1],'Input argument f must have dimensions (m * 1) where m = NUMEL(f).'); end
if m==0 && numel(f)>0, error('Input argument f must have dimensions (m * 1) where m = NUMEL(f).'); end
checkdim(df0dx,[n 1],'Input argument df0dx must have dimensions (n * 1) where n = NUMEL(x).');
if m>0, checkdim(dfdx,[m n],'Input argument dfdx must have dimensions (m * n) where m = NUMEL(f) and n = NUMEL(x).'); end
if m==0 && numel(dfdx)>0, error('Input argument dfdx must have dimensions (m * n) where m = NUMEL(f) and n = NUMEL(x).'); end
checkdim(a0,[1 1],'Input argument a0 must have dimensions (1 * 1).');
checkdim(a,[m 1],'Input argument a must have dimensions (1 * 1) or (m * 1) where m = NUMEL(f).');
checkdim(c,[m 1],'Input argument c must have dimensions (1 * 1) or (m * 1) where m = NUMEL(f).');
checkdim(d,[m 1],'Input argument d must have dimensions (1 * 1) or (m * 1) where m = NUMEL(f).');
checkdim(move,[n 1],'Input argument move must have dimensions (1 * 1) or (n * 1) where n = NUMEL(x).');
checkdim(low,[n 1],'Input argument low must have dimensions (0 * 0), (1 * 1), or (n * 1) where n = NUMEL(x).'); 
checkdim(upp,[n 1],'Input argument upp must have dimensions (0 * 0), (1 * 1), or (n * 1) where n = NUMEL(x).'); 
checkdim(sinit,[n 1],'Input argument sinit must have dimensions (1 * 1) or (n * 1) where n = NUMEL(x).');
checkdim(sincr,[n 1],'Input argument sincr must have dimensions (1 * 1) or (n * 1) where n = NUMEL(x).');
checkdim(sdecr,[n 1],'Input argument sdecr must have dimensions (1 * 1) or (n * 1) where n = NUMEL(x).');
checkdim(mu,[n 1],'Input argument mu must have dimensions (1 * 1) or (n * 1) where n = NUMEL(x).');

% INITIALIZE INTERNAL MMA PARAMETERS ON FIRST ITERATION
if isempty(mmaparams)
  mmaparams = struct;
  mmaparams.iter = 0;
  mmaparams.xold1 = nan(n,1);
  mmaparams.xold2 = nan(n,1);
  mmaparams.low = nan(n,1);
  mmaparams.upp = nan(n,1);
  mmaparams.history.time = [];
  mmaparams.history.iter = [];
  mmaparams.history.funevals = [];
  mmaparams.history.x = [];
  mmaparams.history.f0 = [];
  mmaparams.history.f = [];
end

% PREPARE INPUT PARAMETERS FOR SVANBERG'S MMASUB FUNCTION
iter = mmaparams.iter;
xold1 = mmaparams.xold1;
xold2 = mmaparams.xold2;
low0 = low;
upp0 = upp;
low = mmaparams.low;
upp = mmaparams.upp;
history = mmaparams.history;

% PERFORM MMA ITERATION
iter = iter+1;
[xnew,y,z,lambda,xi,eta,mu,zeta,s,low,upp,alpha,beta,p0,q0,r0,p,q,r] = ...
mmasub(m,n,iter,x,xmin,xmax,xold1,xold2,f0,df0dx,f,dfdx,low,upp,a0,a,c,d,move,low0,upp0,sinit,sincr,sdecr,mu);
xold2 = xold1;
xold1 = x;
change = max(abs(xold1-xold2));

% PRINT CONVERGENCE HISTORY INFORMATION
time = now;
if ~silent
  fprintf('[%s]  iter: %5i | change: % 11.5f | f0: % 11.5f | fmax: % 11.5f\n',datestr(time,'HH:MM:SS'),iter,change,f0,max(f));
end

% COLLECT LAGRANGE MULTIPLIERS AND SLACK VARIABLES
lmult = struct;
lmult.lambda = lambda;
lmult.xi = xi;
lmult.eta = eta;
lmult.mu = mu;
lmult.zeta = zeta;
lmult.s = s;

% COLLECT CONVERGENCE HISTORY INFORMATION
if nargout>=8
  history.time(:,iter) = time;
  history.iter(:,iter) = iter;
  history.funevals(:,iter) = iter;
  if n<=100, history.x(:,iter) = x; end
  history.change(:,iter) = change;
  history.f0(:,iter) = f0;
  if m<=100, history.f(:,iter) = f; else history.f(:,iter) = max(f); end
end

% RETURN CHANGE=MAX(XMAX-XMIN) INSTEAD OF CHANGE=NAN FOR THE FIRST ITERATION SUCH THAT THE STOP CRITERION IS NOT TRIGGERED YET
if iter==1, change=max(xmax-xmin); end

% COLLECT INTERNAL MMA PARAMETERS
mmaparams = struct;
mmaparams.iter = iter;
mmaparams.xold2 = xold2;
mmaparams.xold1 = xold1;
mmaparams.low = low;
mmaparams.upp = upp;
mmaparams.history = history;

% COLLECT MMA SUBPROBLEM INFORMATION
subp = struct;
subp.p0 = p0;
subp.q0 = q0;
subp.r0 = r0;
subp.p = p;
subp.q = q;
subp.r = r;
subp.low = low;
subp.upp = upp;
subp.low(~isnan(low0)) = low0(~isnan(low0));
subp.upp(~isnan(upp0)) = upp0(~isnan(upp0));
subp.alpha = alpha;
subp.beta = beta;


% SUBFUNCTION CHECKDIM - CHECK INPUT ARGUMENT DIMENSIONS
function checkdim(x,dim,err)
if length(size(x))~=length(dim) || any(size(x)~=dim), error(err); end

%-------------------------------------------------------------
%
%    Copyright (C) 2007, 2008 Krister Svanberg
%
%    This file, mmasub.m, is part of GCMMA-MMA-code.
%
%    GCMMA-MMA-code is free software; you can redistribute it and/or
%    modify it under the terms of the GNU General Public License as
%    published by the Free Software Foundation; either version 3 of
%    the License, or (at your option) any later version.
%
%    This code is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    (file COPYING) along with this file.  If not, see
%    <http://www.gnu.org/licenses/>.
%
%    You should have received a file README along with this file,
%    containing contact information.  If not, see
%    <http://www.smoptit.se/> or e-mail mmainfo@smoptit.se or krille@math.kth.se.
%
%    Version September 2007 (and a small change August 2008)
%
%
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp,alfa,beta,p0,q0,r0,P,Q,r] = ...     % MS: added alfa,beta,p0,q0,r0,P,Q,r as output arguments
mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d,move,low0,upp0,asyinit,asyincr,asydecr,albefa);    % MS: added move,low0,upp0,asyinit,asyincr,asydecr,albefa as input arguments
%
%    This function mmasub performs one MMA-iteration, aimed at
%    solving the nonlinear programming problem:
%
%      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
%    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
%                xmin_j <= x_j <= xmax_j,    j = 1,...,n
%                z >= 0,   y_i >= 0,         i = 1,...,m
%*** INPUT:
%
%   m    = The number of general constraints.
%   n    = The number of variables x_j.
%  iter  = Current iteration number ( =1 the first time mmasub is called).
%  xval  = Column vector with the current values of the variables x_j.
%  xmin  = Column vector with the lower bounds for the variables x_j.
%  xmax  = Column vector with the upper bounds for the variables x_j.
%  xold1 = xval, one iteration ago (provided that iter>1).
%  xold2 = xval, two iterations ago (provided that iter>2).
%  f0val = The value of the objective function f_0 at xval.
%  df0dx = Column vector with the derivatives of the objective function
%          f_0 with respect to the variables x_j, calculated at xval.
%  fval  = Column vector with the values of the constraint functions f_i,
%          calculated at xval.
%  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
%          f_i with respect to the variables x_j, calculated at xval.
%          dfdx(i,j) = the derivative of f_i with respect to x_j.
%  low   = Column vector with the lower asymptotes from the previous
%          iteration (provided that iter>1).
%  upp   = Column vector with the upper asymptotes from the previous
%          iteration (provided that iter>1).
%  a0    = The constants a_0 in the term a_0*z.
%  a     = Column vector with the constants a_i in the terms a_i*z.
%  c     = Column vector with the constants c_i in the terms c_i*y_i.
%  d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
%
%*** OUTPUT:
%
%  xmma  = Column vector with the optimal values of the variables x_j
%          in the current MMA subproblem.
%  ymma  = Column vector with the optimal values of the variables y_i
%          in the current MMA subproblem.
%  zmma  = Scalar with the optimal value of the variable z
%          in the current MMA subproblem.
%  lam   = Lagrange multipliers for the m general MMA constraints.
%  xsi   = Lagrange multipliers for the n constraints alfa_j - x_j <= 0.
%  eta   = Lagrange multipliers for the n constraints x_j - beta_j <= 0.
%   mu   = Lagrange multipliers for the m constraints -y_i <= 0.
%  zet   = Lagrange multiplier for the single constraint -z <= 0.
%   s    = Slack variables for the m general MMA constraints.
%  low   = Column vector with the lower asymptotes, calculated and used
%          in the current MMA subproblem.
%  upp   = Column vector with the upper asymptotes, calculated and used
%          in the current MMA subproblem.
%
%epsimin = sqrt(m+n)*10^(-9);
epsimin = 10^(-7);
raa0 = 0.00001;
% move = 0.5;       % MS: parameter is now an input argument
% albefa = 0.1;     % MS: parameter is now an input argument
% asyinit = 0.5;    % MS: parameter is now an input argument
% asyincr = 1.2;    % MS: parameter is now an input argument
% asydecr = 0.7;    % MS: parameter is now an input argument
eeen = ones(n,1);
eeem = ones(m,1);
zeron = zeros(n,1);
albefa(isinf(low0)&isinf(upp0)) = 0; % MS: do not use mu/albefa in the case of SLP
albefa((low0==0)&isinf(upp0)) = 0; % MS: do not use mu/albefa in the case of CONLIN

% Calculation of the asymptotes low and upp :
if iter < 2.5
  low = xval - asyinit.*(xmax-xmin); % MS: .* instead of * as asyinit now is a vector
  upp = xval + asyinit.*(xmax-xmin); % MS: .* instead of * as asyinit now is a vector
else
  zzz = (xval-xold1).*(xold1-xold2);
  factor = eeen;
  factor(find(zzz > 0)) = asyincr(find(zzz > 0)); % MS: added indices (find(...)) at RHS as asyincr now is a vector
  factor(find(zzz < 0)) = asydecr(find(zzz < 0)); % MS: added indices (find(...)) at RHS as asydecr now is a vector
  low = xval - factor.*(xold1 - low);
  upp = xval + factor.*(upp - xold1);
  lowmin = xval - 10*(xmax-xmin);
  lowmax = xval - 0.01*(xmax-xmin);
  uppmin = xval + 0.01*(xmax-xmin);
  uppmax = xval + 10*(xmax-xmin);
  low = max(low,lowmin);
  low = min(low,lowmax);
  upp = min(upp,uppmax);
  upp = max(upp,uppmin);
end

% MS: do not allow asymptotes to move outside box constraints initially
if iter < 2.5
  low = max(xmin,low);
  upp = min(xmax,upp);
end

% MS: if user specified values low0 and upp0 are passed, use these instead
low(~isnan(low0)) = low0(~isnan(low0));
upp(~isnan(upp0)) = upp0(~isnan(upp0));

% Calculation of the bounds alfa and beta :

zzz1 = low + albefa.*(xval-low); % MS: .* instead of * as albefa now is a vector
zzz2 = xval - move.*(xmax-xmin); % MS: .* instead of * as move now is a vector
zzz  = max(zzz1,zzz2);
alfa = max(zzz,xmin);
zzz1 = upp - albefa.*(upp-xval); % MS: .* instead of * as albefa now is a vector
zzz2 = xval + move.*(xmax-xmin); % MS: .* instead of * as move now is a vector
zzz  = min(zzz1,zzz2);
beta = min(zzz,xmax);

% MS: replace INF asymptotes with large numbers
low(isinf(low)) = sign(low(isinf(low))).*(beta(isinf(low))-alfa(isinf(low)))*1e4;
upp(isinf(upp)) = sign(upp(isinf(upp))).*(beta(isinf(upp))-alfa(isinf(upp)))*1e4;

% Calculations of p0, q0, P, Q and b.

xmami = xmax-xmin;
xmamieps = 0.00001*eeen;
xmami = max(xmami,xmamieps);
xmamiinv = eeen./xmami;
ux1 = upp-xval;
ux2 = ux1.*ux1;
xl1 = xval-low;
xl2 = xl1.*xl1;
uxinv = eeen./ux1;
xlinv = eeen./xl1;
%
p0 = zeron;
q0 = zeron;
p0 = max(df0dx,0);
q0 = max(-df0dx,0);
pq0 = 0.001*(p0 + q0) + raa0*xmamiinv;
pq0(isinf(low0)|isinf(upp0)) = 0; % MS: do not add additional terms from 2007 in the case of SLP or CONLIN
p0 = p0 + pq0;
q0 = q0 + pq0;
p0 = p0.*ux2;
q0 = q0.*xl2;
%
P = sparse(m,n);
Q = sparse(m,n);
P = max(dfdx,0);
Q = max(-dfdx,0);
PQ = 0.001*(P + Q) + raa0*eeem*xmamiinv';
PQ(:,isinf(low0)|isinf(upp0)) = 0; % MS: do not add additional terms from 2007 in the case of SLP or CONLIN
P = P + PQ;
Q = Q + PQ;
P = P * spdiags(ux2,0,n,n);
Q = Q * spdiags(xl2,0,n,n);
b = P*uxinv + Q*xlinv - fval ;
%
% MS: compute constant terms r0 and r in the convex approximation
b0 = p0.'*uxinv + q0.'*xlinv - f0val;
r0 = -b0;
r = -b;

%%% Solving the subproblem by a primal-dual Newton method
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);


%-------------------------------------------------------------
%
%    Copyright (C) 2006 Krister Svanberg
%
%    This file, subsolv.m, is part of GCMMA-MMA-code.
%
%    GCMMA-MMA-code is free software; you can redistribute it and/or
%    modify it under the terms of the GNU General Public License as
%    published by the Free Software Foundation; either version 3 of
%    the License, or (at your option) any later version.
%
%    This code is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    (file COPYING) along with this file.  If not, see
%    <http://www.gnu.org/licenses/>.
%
%    You should have received a file README along with this file,
%    containing contact information.  If not, see
%    <http://www.smoptit.se/> or e-mail mmainfo@smoptit.se or krille@math.kth.se.
%
%    Version Dec 2006.
%
%
function [xmma,ymma,zmma,lamma,xsimma,etamma,mumma,zetmma,smma] = ...
subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);
%
% This function subsolv solves the MMA subproblem:
%
% minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
%          + SUM[ ci*yi + 0.5*di*(yi)^2 ],
%
% subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
%            alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
%
% Input:  m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
% Output: xmma,ymma,zmma, slack variables and Lagrange multiplers.
%

% MS: added code for unconstrained problems
if m==0
  f0alfa = sum(p0./(upp-alfa)+q0./(alfa-low));
  f0beta = sum(p0./(upp-beta)+q0./(beta-low));
  kalfa = f0alfa<=f0beta;
  kbeta = f0beta<f0alfa;
  xmma = zeros(n,1);
  xmma(kalfa) = alfa(kalfa);
  xmma(kbeta) = beta(kbeta);
  df0dx = p0./(upp-xmma).^2-q0./(xmma-low).^2;
  xsimma = zeros(n,1);
  xsimma(kalfa) = df0dx(kalfa);
  etamma = zeros(n,1);
  etamma(kbeta) = -df0dx(kbeta);
  ymma = [];
  if a0>=0,
    zmma=0;
    zetmma = a0;
  else
    zmma=inf;
    zetmma = 0;
  end
  lamma = [];
  mumma = [];
  smma = [];
  return;
end

een = ones(n,1);
eem = ones(m,1);
epsi = 1;
epsvecn = epsi*een;
epsvecm = epsi*eem;
x = 0.5*(alfa+beta);
y = eem;
z = 1;
lam = eem;
xsi = een./(x-alfa);
xsi = max(xsi,een);
eta = een./(beta-x);
eta = max(eta,een);
mu  = max(eem,0.5*c);
zet = 1;
s = eem;
itera = 0;

while epsi > epsimin
  epsvecn = epsi*een;
  epsvecm = epsi*eem;
  ux1 = upp-x;
  xl1 = x-low;
  ux2 = ux1.*ux1;
  xl2 = xl1.*xl1;
  uxinv1 = een./ux1;
  xlinv1 = een./xl1;
  plam = p0 + P'*lam ;
  qlam = q0 + Q'*lam ;
  gvec = P*uxinv1 + Q*xlinv1;
  dpsidx = plam./ux2 - qlam./xl2 ;
  rex = dpsidx - xsi + eta;
  rey = c + d.*y - mu - lam;
  rez = a0 - zet - a'*lam;
  relam = gvec - a*z - y + s - b;
  rexsi = xsi.*(x-alfa) - epsvecn;
  reeta = eta.*(beta-x) - epsvecn;
  remu = mu.*y - epsvecm;
  rezet = zet*z - epsi;
  res = lam.*s - epsvecm;
  residu1 = [rex' rey' rez]';
  residu2 = [relam' rexsi' reeta' remu' rezet res']';
  residu = [residu1' residu2']';
  residunorm = sqrt(residu'*residu);
  residumax = max(abs(residu));
  ittt = 0;
  while residumax > 0.9*epsi & ittt < 200
    ittt=ittt + 1;
    itera=itera + 1;
    ux1 = upp-x;
    xl1 = x-low;
    ux2 = ux1.*ux1;
    xl2 = xl1.*xl1;
    ux3 = ux1.*ux2;
    xl3 = xl1.*xl2;
    uxinv1 = een./ux1;
    xlinv1 = een./xl1;
    uxinv2 = een./ux2;
    xlinv2 = een./xl2;
    plam = p0 + P'*lam ;
    qlam = q0 + Q'*lam ;
    gvec = P*uxinv1 + Q*xlinv1;
    GG = P*(uxinv2.*speye(n,n)) - Q*(xlinv2.*speye(n,n));
    dpsidx = plam./ux2 - qlam./xl2 ;
    delx = dpsidx - epsvecn./(x-alfa) + epsvecn./(beta-x);
    dely = c + d.*y - lam - epsvecm./y;
    delz = a0 - a'*lam - epsi/z;
    dellam = gvec - a*z - y - b + epsvecm./lam;
    diagx = plam./ux3 + qlam./xl3;
    diagx = 2*diagx + xsi./(x-alfa) + eta./(beta-x);
    diagxinv = een./diagx;
    diagy = d + mu./y;
    diagyinv = eem./diagy;
    diaglam = s./lam;
    diaglamyi = diaglam+diagyinv;
    if m < n
      blam = dellam + dely./diagy - GG*(delx./diagx);
      bb = [blam' delz]';
      Alam = (diaglamyi.*speye(m,m)) + GG*((diagxinv.*speye(n,n))*GG');
      AA = [Alam     a
            a'    -zet/z ];
      solut = AA\bb;
      dlam = solut(1:m);
      dz = solut(m+1);
      dx = -delx./diagx - (GG'*dlam)./diagx;
    else
      diaglamyiinv = eem./diaglamyi;
      dellamyi = dellam + dely./diagy;
      Axx = (diagx.*speye(n,n)) + GG'*((diaglamyiinv.*speye(m,m))*GG);
      azz = zet/z + a'*(a./diaglamyi);
      axz = -GG'*(a./diaglamyi);
      bx = delx + GG'*(dellamyi./diaglamyi);
      bz  = delz - a'*(dellamyi./diaglamyi);
      AA = [Axx   axz
            axz'  azz ];
      bb = [-bx' -bz]';
      solut = AA\bb;
      dx  = solut(1:n);
      dz = solut(n+1);
      dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
    end
%
    dy = -dely./diagy + dlam./diagy;
    dxsi = -xsi + epsvecn./(x-alfa) - (xsi.*dx)./(x-alfa);
    deta = -eta + epsvecn./(beta-x) + (eta.*dx)./(beta-x);
    dmu  = -mu + epsvecm./y - (mu.*dy)./y;
    dzet = -zet + epsi/z - zet*dz/z;
    ds   = -s + epsvecm./lam - (s.*dlam)./lam;
    xx  = [ y'  z  lam'  xsi'  eta'  mu'  zet  s']';
    dxx = [dy' dz dlam' dxsi' deta' dmu' dzet ds']';
%
    stepxx = -1.01*dxx./xx;
    stmxx  = max(stepxx);
    stepalfa = -1.01*dx./(x-alfa);
    stmalfa = max(stepalfa);
    stepbeta = 1.01*dx./(beta-x);
    stmbeta = max(stepbeta);
    stmalbe  = max(stmalfa,stmbeta);
    stmalbexx = max(stmalbe,stmxx);
    stminv = max(stmalbexx,1);
    steg = 1/stminv;
%
    xold   =   x;
    yold   =   y;
    zold   =   z;
    lamold =  lam;
    xsiold =  xsi;
    etaold =  eta;
    muold  =  mu;
    zetold =  zet;
    sold   =   s;
%
    itto = 0;
    resinew = 2*residunorm;
    while resinew > residunorm & itto < 50
    itto = itto+1;
    x   =   xold + steg*dx;
    y   =   yold + steg*dy;
    z   =   zold + steg*dz;
    lam = lamold + steg*dlam;
    xsi = xsiold + steg*dxsi;
    eta = etaold + steg*deta;
    mu  = muold  + steg*dmu;
    zet = zetold + steg*dzet;
    s   =   sold + steg*ds;
    ux1 = upp-x;
    xl1 = x-low;
    ux2 = ux1.*ux1;
    xl2 = xl1.*xl1;
    uxinv1 = een./ux1;
    xlinv1 = een./xl1;
    plam = p0 + P'*lam ;
    qlam = q0 + Q'*lam ;
    gvec = P*uxinv1 + Q*xlinv1;
    dpsidx = plam./ux2 - qlam./xl2 ;
    rex = dpsidx - xsi + eta;
    rey = c + d.*y - mu - lam;
    rez = a0 - zet - a'*lam;
    relam = gvec - a*z - y + s - b;
    rexsi = xsi.*(x-alfa) - epsvecn;
    reeta = eta.*(beta-x) - epsvecn;
    remu = mu.*y - epsvecm;
    rezet = zet*z - epsi;
    res = lam.*s - epsvecm;
    residu1 = [rex' rey' rez]';
    residu2 = [relam' rexsi' reeta' remu' rezet res']';
    residu = [residu1' residu2']';
    resinew = sqrt(residu'*residu);
    steg = steg/2;
    end
  residunorm=resinew;
  residumax = max(abs(residu));
  steg = 2*steg;
  end
  % if ittt > 198  % MS: do not show warning
  %   epsi
  %   ittt
  % end

epsi = 0.1*epsi;
end
xmma   =   x;
ymma   =   y;
zmma   =   z;
lamma =  lam;
xsimma =  xsi;
etamma =  eta;
mumma  =  mu;
zetmma =  zet;
smma   =   s;
%-------------------------------------------------------------
