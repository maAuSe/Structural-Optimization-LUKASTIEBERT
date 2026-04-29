function [xnew,y,z,lmult,mmaparams,subp,change,history] = gcmma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams,move,low,upp,sinit,sincr,sdecr,mu,a0,a,c,d,arg1,arg2)

%GCMMA   Globally convergent method of moving asymptotes.
%
%   [xnew,y,z,lmult,mmaparams,subp,change,history] = ...
%   GCMMA(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams,move,low,upp,sinit,sincr,sdecr,mu,a0,a,c,d)
%   performs one GCMMA-iteration, aimed at solving the nonlinear programming problem:
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
%                     Not needed for inner iterations - in this case pass any value.
%   dfdx              Derivatives of the constraints at x (m * n).
%                     Not needed for inner iterations - in this case pass any value.
%   mmaparams         Internal parameters controlling the moving asymptotes and the
%                     level of conservatism in the convex approximations.  Use []
%                     in the first iteration and return the corresponding output
%                     argument in the following iterations.
%   move              Move limit (n * 1) or (1 * 1). Default: 0.2.
%   low               User-specified values for the lower asymptotes (n * 1) or (1 * 1). Default: [].
%   upp               User-specified values for the upper asymptotes (n * 1) or (1 * 1). Default: [].
%                     Specify [] or nan to compute asymptotes automatically.
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
%   xnew              Optimal values of the variables x in the current GCMMA subproblem (n * 1).
%   y                 Optimal values of the variables y in the current GCMMA subproblem (n * 1).
%   z                 Optimal values of the variables z in the current GCMMA subproblem (n * 1).
%   lmult.lambda      Lagrange multipliers for the general GCMMA constraints (m * 1).
%   lmult.xi          Lagrange multipliers for the constraints alfa_j - x_j <= 0 (n * 1).
%   lmult.eta         Lagrange multipliers for the constraints x_j - beta_j <= 0 (n * 1).
%   lmult.mu          Lagrange multipliers for the constraints -y_i <= 0 (m * 1).
%   lmult.zeta        Lagrange multiplier for the constraint -z <= 0 (1 * 1).
%   lmult.s           Slack variables for the general GCMMA constraints (m * 1).
%   mmaparams         Updated values of the internal GCMMA parameters.
%   subp.p0           Coefficients of the upper asymptote terms in the objective function approximation (n * 1).
%   subp.q0           Coefficients of the lower asymptote terms in the objective function approximation (n * 1).
%   subp.r0           Constant term in the objective function approximation (1 * 1).
%   subp.p            Coefficients of the upper asymptote terms in the constraint approximations (m * n).
%   subp.q            Coefficients of the lower asymptote terms in the constraint approximations (m * n).
%   subp.r            Constant term in the constraint approximations (m * 1).
%   subp.low          Values of the lower asymptotes (n * 1).
%   subp.upp          Values of the upper asymptotes (n * 1).
%   subp.alpha        Lower bounds for x used in the current GCMMA subproblem (n * 1).
%   subp.beta         Upper bounds for x used in the current GCMMA subproblem (n * 1).
%   subp.f0           Value of the objective function approximation at xnew (1 * 1).
%   subp.f            Values of the constraint function approximations at xnew (1 * 1).
%   subp.f0e          Equal to subp.f0+1e-7, used to determine if an inner iteration is needed next.
%   subp.fe           Equal to subp.f+1e-7, used to determine if an inner iteration is needed next.
%                     An inner iteration is needed next if the approximation is not conservative at xnew.
%                     The approximation is conservative if f0<=subp.f0e && all(f<=subp.fe)
%   change            Design change, defined as MAX(ABS(x-xold)) if the current iteration is an outer loop
%                     iteration, where x is the design passed as input in the current iteration, and xold is
%                     the design passed as input in the previous outer loop iteration. If the current
%                     iteration is an inner loop iteration, INF is returned instead.  This indicates that
%                     the previous iteration was not conservative, and that more iterations are needed.
%   history.time      Time at outer iteration i (1 * N).
%   history.outeriter Outer iteration counter (1 * N).
%   history.inneriter Number of inner iterations performed between outer iteration i-1 and i (1 * N).
%   history.funevals  Total number of function evaluations performed at outer iteration i (1 * N).
%   history.x         Values of variables x at outer iteration i (n * N).
%                     If the number of design variables exceeds 100, an empty matrix is returned.
%   history.change    Design change, defined as above (1 * N).
%   history.f0        Value of objective function at outer iteration i (1 * N).
%   history.f         Value of constraints at outer iteration i (m * N).
%                     If the number of constraints exceeds 100, MAX(f) is returned instead (1 * N).
%
%   The following information is printed as output to the screen:
%
%   outeriter:        Outer iteration counter.
%   inneriter:        Number of inner iterations performed after the previous outer iteration.
%   funevals:         Total number of function evaluations performed (sum of outer and inner iterations).
%   change:           Design change, defined as above.
%   f0:               Value of the objective function at x.
%   fmax:             Value of the most critical constraint at x, i.e. MAX(f).
%
%   This function follows Svanberg's implementation (www.smoptit.se) very minor modifications.
%   It includes an optional adaptive update of the initial conservativity parameter (raa) based on the 
%   number of inner iterations in the previous outer step.  This aims to reduce repeated inner iterations.  
%   Activate with GCMMA(...,'innercontrol').
%
%   Use GCMMA(...,'silent') to avoid printing output.

% Mattias Schevenels
% January 2026

% Modifications with respect to Svanberg's code are denoted by the tag MS.

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
if exist('arg2','var') && strcmpi(arg2,'silent'), silent=true; clear('arg2'); end

% DETERMINE INNERCONTROL STATUS
innercontrol = false;
if exist('move','var') && strcmpi(move,'innercontrol'), innercontrol=true; clear('move'); end
if exist('low','var') && strcmpi(low,'innercontrol'), innercontrol=true; clear('low'); end
if exist('upp','var') && strcmpi(upp,'innercontrol'), innercontrol=true; clear('upp'); end
if exist('sinit','var') && strcmpi(sinit,'innercontrol'), innercontrol=true; clear('sinit'); end
if exist('sincr','var') && strcmpi(sincr,'innercontrol'), innercontrol=true; clear('sincr'); end
if exist('sdecr','var') && strcmpi(sdecr,'innercontrol'), innercontrol=true; clear('sdecr'); end
if exist('mu','var') && strcmpi(mu,'innercontrol'), innercontrol=true; clear('mu'); end
if exist('a0','var') && strcmpi(a0,'innercontrol'), innercontrol=true; clear('a0'); end
if exist('a','var') && strcmpi(a,'innercontrol'), innercontrol=true; clear('a'); end
if exist('c','var') && strcmpi(c,'innercontrol'), innercontrol=true; clear('c'); end
if exist('d','var') && strcmpi(d,'innercontrol'), innercontrol=true; clear('d'); end
if exist('arg1','var') && strcmpi(arg1,'innercontrol'), innercontrol=true; clear('arg1'); end
if exist('arg2','var') && strcmpi(arg2,'innercontrol'), innercontrol=true; clear('arg2'); end

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

% CHECK IF USER ATTEMPTS TO USE CONLIN/SLP
if any(isinf(low)) || any(isinf(upp)) error('Infinite asymptotes are not supported in GCMMA; use MMA for SLP/CONLIN.'); end

% FIXED GCMMA PARAMETERS
raa0eps = 1e-6;
raaeps = repmat(1e-6,m,1);
epsimin = 1e-7;

% INITIALIZE INTERNAL MMA PARAMETERS ON FIRST ITERATION
if isempty(mmaparams)
  mmaparams = struct;
  mmaparams.iter = 0;
  mmaparams.outeriter = 0;
  mmaparams.inneriter = 0;
  mmaparams.xold1 = nan(n,1);
  mmaparams.xold2 = nan(n,1);
  mmaparams.low = nan(n,1);
  mmaparams.upp = nan(n,1);
  mmaparams.f0app = inf;
  mmaparams.fapp = inf(m,1);
  mmaparams.raa0 = nan;
  mmaparams.raa = nan(m,1);
  mmaparams.raacoef = 0.1;
  mmaparams.f0old = nan;
  mmaparams.df0dxold = nan(n,1);
  mmaparams.fold = nan(m,1);
  mmaparams.dfdxold = nan(m,n);
  mmaparams.history.time = [];
  mmaparams.history.outeriter = [];
  mmaparams.history.inneriter = [];
  mmaparams.history.funevals = [];
  mmaparams.history.x = [];
  mmaparams.history.f0 = [];
  mmaparams.history.f = [];
end

% READ INTERNAL MMA PARAMETERS
iter = mmaparams.iter;
outeriter = mmaparams.outeriter;
inneriter = mmaparams.inneriter;
xold1 = mmaparams.xold1;
xold2 = mmaparams.xold2;
low0 = low;
upp0 = upp;
low = mmaparams.low;
upp = mmaparams.upp;
f0app = mmaparams.f0app;
fapp = mmaparams.fapp;
f0old = mmaparams.f0old;
df0dxold = mmaparams.df0dxold;
fold = mmaparams.fold;
dfdxold = mmaparams.dfdxold;
raa0 = mmaparams.raa0;
raa = mmaparams.raa;
raacoef = mmaparams.raacoef;
history = mmaparams.history;

% INCREASE ITERATION COUNTER
iter = iter+1;

% DETERMINE WHETHER PREVIOUS ITERATION WAS CONSERVATIVE
if f0<=f0app+epsimin & all(f<=fapp+epsimin);

  % PERFORM OUTER MMA ITERATION
  outeriter = outeriter+1;
  [low,upp,raa0,raa,raacoef] = asymp(outeriter,n,x,xold1,xold2,xmin,xmax,low,upp,[],[],raa0eps,raaeps,df0dx,dfdx,low0,upp0,sinit,sincr,sdecr,raacoef,innercontrol,inneriter);
  [xnew,y,z,lambda,xi,eta,mu,zeta,s,f0app,fapp,alpha,beta,p0,q0,r0,p,q,r] = gcmmasub(m,n,outeriter,epsimin,x,xmin,xmax,low,upp,raa0,raa,f0,df0dx,f,dfdx,a0,a,c,d,move,low0,upp0,mu);
  xold2 = xold1;
  xold1 = x;
  f0old = f0;
  fold = f;
  df0dxold = df0dx;
  dfdxold = dfdx;
  change = max(abs(xold1-xold2));

  % PRINT CONVERGENCE HISTORY INFORMATION
  time = now;
  if ~silent
    fprintf('[%s]  outeriter: %5i | inneriter: %5i | funevals: %5i | change: % 11.5f | f0: % 11.5f | fmax: % 11.5f\n',datestr(time,'HH:MM:SS'),outeriter,inneriter,iter,change,f0,max(f));
  end

  % COLLECT CONVERGENCE HISTORY INFORMATION
  if nargout>=8
    history.time(:,outeriter) = time;
    history.outeriter(:,outeriter) = outeriter;
    history.inneriter(:,outeriter) = inneriter;
    history.funevals(:,outeriter) = iter;
    if n<=100, history.x(:,outeriter) = x; end
    history.change(:,outeriter) = change;
    history.f0(:,outeriter) = f0;
    if m<=100, history.f(:,outeriter) = f; else history.f(:,outeriter) = max(f); end
  end

  % RESET INNER ITERATION COUNTER
  inneriter = 0;
else

  % PERFORM INNER MMA ITERATION
  inneriter = inneriter+1;
  [raa0,raa] = raaupdate(x,xold1,xmin,xmax,low,upp,f0,f,f0app,fapp,raa0,raa,raa0eps,raaeps,epsimin);
  [xnew,y,z,lambda,xi,eta,mu,zeta,s,f0app,fapp,alpha,beta,p0,q0,r0,p,q,r] = gcmmasub(m,n,outeriter,epsimin,xold1,xmin,xmax,low,upp,raa0,raa,f0old,df0dxold,fold,dfdxold,a0,a,c,d,move,low0,upp0,mu);
  change = inf;

end

% RETURN CHANGE=MAX(XMAX-XMIN) INSTEAD OF CHANGE=NAN FOR THE FIRST ITERATION SUCH THAT THE STOP CRITERION IS NOT TRIGGERED YET
if iter==1, change=max(xmax-xmin); end

% LIMIT THE NUMBER OF CONSECUTIVE INNER ITERATIONS TO 16
if inneriter>=16
  f0app = inf;
  fapp = inf(m,1);
end


% COLLECT LAGRANGE MULTIPLIERS AND SLACK VARIABLES
lmult = struct;
lmult.lambda = lambda;
lmult.xi = xi;
lmult.eta = eta;
lmult.mu = mu;
lmult.zeta = zeta;
lmult.s = s;

% COLLECT INTERNAL MMA PARAMETERS
mmaparams = struct;
mmaparams.iter = iter;
mmaparams.outeriter = outeriter;
mmaparams.inneriter = inneriter;
mmaparams.xold2 = xold2;
mmaparams.xold1 = xold1;
mmaparams.low = low;
mmaparams.upp = upp;
mmaparams.f0app = f0app;
mmaparams.fapp = fapp;
mmaparams.f0old = f0old;
mmaparams.df0dxold = df0dxold;
mmaparams.fold = fold;
mmaparams.dfdxold = dfdxold;
mmaparams.raa0 = raa0;
mmaparams.raa = raa;
mmaparams.raacoef = raacoef;
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
subp.f0 = f0app;
subp.f = fapp;
subp.f0e = f0app+epsimin;
subp.fe = fapp+epsimin;

% SUBFUNCTION CHECKDIM - CHECK INPUT ARGUMENT DIMENSIONS
function checkdim(x,dim,err)
if length(size(x))~=length(dim) || any(size(x)~=dim), error(err); end

%-------------------------------------------------------------
%
%    Copyright (C) 2007 Krister Svanberg
%
%    This file, asymp.m, is part of GCMMA-MMA-code.
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
%------
%
%  Values on the parameters raa0, raa, low and upp are
%  calculated in the beginning of each outer iteration.
%
function [low,upp,raa0,raa,raacoef] = ...
asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp, ...
      raa0,raa,raa0eps,raaeps,df0dx,dfdx,low0,upp0,asyinit,asyincr,asydecr,raacoef,innercontrol,inneriter); % MS: added low0,upp0,asyinit,asyincr,asydecr,raacoef,innercontrol,inneriter as input arguments
%
eeen=ones(n,1);
xmami = xmax - xmin;
xmamieps = 0.00001*eeen;
xmami = max(xmami,xmamieps);
raa0 = abs(df0dx)'*xmami;
%raa0 = max(raa0eps,(0.1/n)*raa0); % MS: replaced with if-statement below 
raa  = abs(dfdx)*xmami;
%raa  = max(raaeps,(0.1/n)*raa); % MS: replaced with if-statement below

% MS: ADAPTIVE RULE FOR INITIAL RAA/RAA0 TO CONTROL THE NUMBER OF INNER ITERATIONS
if innercontrol                           % MS
  if inneriter>0                          % MS
    raacoef = min(100,2*raacoef);         % MS
  else                                    % MS
    raacoef = max(0.1,0.5^(1/5)*raacoef); % MS (5 of these steps needed to compensate for increase of conservativeness)
  end                                     % MS
  raa0 = max(raa0eps,(raacoef/n)*raa0);   % MS
  raa  = max(raaeps,(raacoef/n)*raa);     % MS
else                                      % standard GCMMA
  raa0 = max(raa0eps,(0.1/n)*raa0);       % standard GCMMA
  raa  = max(raaeps,(0.1/n)*raa);         % standard GCMMA
end                                     

if outeriter < 2.5
  low = xval - asyinit.*xmami; % MS: .* instead of * as asyinit now is a vector 
  upp = xval + asyinit.*xmami; % MS: .* instead of * as asyinit now is a vector
else
  xxx = (xval-xold1).*(xold1-xold2);
  factor = eeen;
  factor(find(xxx > 0)) = asyincr(find(xxx > 0)); % MS: added indices (find(...)) at RHS as asyincr now is a vector
  factor(find(xxx < 0)) = asydecr(find(xxx < 0)); % MS: added indices (find(...)) at RHS as asydecr now is a vector
  low = xval - factor.*(xold1 - low);
  upp = xval + factor.*(upp - xold1);
  lowmin = xval - 10*xmami;
  lowmax = xval - 0.01*xmami;
  uppmin = xval + 0.01*xmami;
  uppmax = xval + 10*xmami;
  low = max(low,lowmin);
  low = min(low,lowmax);
  upp = min(upp,uppmax);
  upp = max(upp,uppmin);
end

% MS: do not allow asymptotes to move outside box constraints initially
if outeriter < 2.5
  low = max(xmin,low);
  upp = min(xmax,upp);
end

% MS: if user specified values low0 and upp0 are passed, use these instead
low(~isnan(low0)) = low0(~isnan(low0));
upp(~isnan(upp0)) = upp0(~isnan(upp0));


%-------------------------------------------------------------
%
%    Copyright (C) 2007 Krister Svanberg
%
%    This file, raaupdate.m, is part of GCMMA-MMA-code.
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
%------
%
%  Version April 2007.
%
%  Values of the parameters raa0 and raa are updated
%  during an inner iteration.
%
function [raa0,raa] = ...
raaupdate(xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew, ...
          f0app,fapp,raa0,raa,raa0eps,raaeps,epsimin);
%
raacofmin = 10^(-12);
eeem = ones(size(raa));
eeen = ones(size(xmma));
xmami = xmax-xmin;
xmamieps = 0.00001*eeen;
xmami = max(xmami,xmamieps);
xxux = (xmma-xval)./(upp-xmma);
xxxl = (xmma-xval)./(xmma-low);
xxul = xxux.*xxxl;
ulxx = (upp-low)./xmami;
raacof = xxul'*ulxx;
raacof = max(raacof,raacofmin);
%
f0appe = f0app + 0.5*epsimin;
if f0valnew > f0appe
  deltaraa0 = (1/raacof)*(f0valnew-f0app);
  zz0 = 1.1*(raa0 + deltaraa0);
  zz0 = min(zz0,10*raa0);
%  zz0 = min(zz0,1000*raa0);
  raa0 = zz0;
end
%
fappe = fapp + 0.5*epsimin*eeem;
fdelta = fvalnew-fappe;
deltaraa = (1/raacof)*(fvalnew-fapp);
zzz = 1.1*(raa + deltaraa);
zzz = min(zzz,10*raa);
%zzz = min(zzz,1000*raa);
raa(find(fdelta > 0)) = zzz(find(fdelta > 0));

%-------------------------------------------------------------
%
%    Copyright (C) 2008 Krister Svanberg
%
%    This file, gcmmasub.m.m, is part of GCMMA-MMA-code.
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
%
%------
%
%  Version Feb 2008.
%
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp,alfa,beta,p0,q0,r0,P,Q,r] = ... % MS: added alfa,beta,p0,q0,r0,P,Q,r as output arguments
gcmmasub(m,n,iter,epsimin,xval,xmin,xmax,low,upp,raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d,move,low0,upp0,albefa); % MS: added move,low0,upp0,albefa as input arguments
%
eeen = ones(n,1);
zeron = zeros(n,1);
%
% Calculations of the bounds alfa and beta.
% albefa = 0.1;  % MS: parameter is now an input argument
% move = 0.5;    % MS: parameter is now an input argument
%
zzz1 = low + albefa.*(xval-low); % MS: .* instead of * as albefa now is a vector
zzz2 = xval - move.*(xmax-xmin); % MS: .* instead of * as move now is a vector
zzz  = max(zzz1,zzz2);
alfa = max(zzz,xmin);
zzz1 = upp - albefa.*(upp-xval); % MS: .* instead of * as albefa now is a vector
zzz2 = xval + move.*(xmax-xmin); % MS: .* instead of * as move now is a vector
zzz  = min(zzz1,zzz2);
beta = min(zzz,xmax);

% Calculations of p0, q0, r0, P, Q, r and b.
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
pq0 = p0 + q0;
p0 = p0 + 0.001*pq0;
q0 = q0 + 0.001*pq0;
p0 = p0 + raa0*xmamiinv;
q0 = q0 + raa0*xmamiinv;
p0 = p0.*ux2;
q0 = q0.*xl2;
r0 = f0val - p0'*uxinv - q0'*xlinv;
%
P = sparse(m,n);
Q = sparse(m,n);
P = max(dfdx,0);
Q = max(-dfdx,0);
PQ = P + Q;
P = P + 0.001*PQ;
Q = Q + 0.001*PQ;
P = P + raa*xmamiinv';
Q = Q + raa*xmamiinv';
P = P * spdiags(ux2,0,n,n);
Q = Q * spdiags(xl2,0,n,n);
r = fval - P*uxinv - Q*xlinv;
b = -r;
%
% Solving the subproblem by a primal-dual Newton method
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);
%
% Calculations of f0app and fapp.
ux1 = upp-xmma;
xl1 = xmma-low;
uxinv = eeen./ux1;
xlinv = eeen./xl1;
f0app = r0 + p0'*uxinv + q0'*xlinv;
fapp  =  r +   P*uxinv +   Q*xlinv;
%

%-------------------------------------------------------------
%    This is the file subsolv.m
%
%    Version Dec 2006.
%    Krister Svanberg <krille@math.kth.se>
%    Department of Mathematics, KTH,
%    SE-10044 Stockholm, Sweden.
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
    GG = P*spdiags(uxinv2,0,n,n) - Q*spdiags(xlinv2,0,n,n);
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
      Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
      AA = [Alam     a
            a'    -zet/z ];
      solut = AA\bb;
      dlam = solut(1:m);
      dz = solut(m+1);
      dx = -delx./diagx - (GG'*dlam)./diagx;
    else
      diaglamyiinv = eem./diaglamyi;
      dellamyi = dellam + dely./diagy;
      Axx = spdiags(diagx,0,n,n) + GG'*spdiags(diaglamyiinv,0,m,m)*GG;
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
