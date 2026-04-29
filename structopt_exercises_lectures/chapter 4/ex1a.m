%EX1A
%
%   Use of SLP, CONLIN, and MMA to solve a 1-D optimization problem.

% Mattias Schevenels
% February 2021

% MMA PARAMETERS
% xmin, xmax:     Box constraints
% x:              Initial value
% move:           Move limit (relative to xmax-xmin)
% low:            Lower asymptote: use -INF for SLP,   0 for CONLIN, [] for MMA
% upp:            Upper asymptote: use  INF for SLP, INF for CONLIN, [] for MMA
% sinit           Asymptotes relative initial position.
% sincr           Factor used to increase the distance between the asymptotes
% sdecr           Factor used to decrease the distance between the asymptotes
% mu              Relative distance to keep away from the asymptotes
preset = 1;
switch preset
  case  1, xmin = 0;   xmax = 10; x =  3; move = inf;  low = -inf; upp = inf; sinit = 0.2;  sincr = 1.2; sdecr = 0.7; mu = 0;  % SLP
  case  2, xmin = 0;   xmax = 10; x = 10; move = inf;  low = -inf; upp = inf; sinit = 0.2;  sincr = 1.2; sdecr = 0.7; mu = 0;  % SLP, different starting point
  case  3, xmin = 0.1; xmax = 10; x = 10; move = inf;  low = -inf; upp = inf; sinit = 0.2;  sincr = 1.2; sdecr = 0.7; mu = 0;  % SLP, increase lower bound
  case  4, xmin = 0;   xmax = 10; x = 10; move = 0.2;  low = -inf; upp = inf; sinit = 0.2;  sincr = 1.2; sdecr = 0.7; mu = 0;  % SLP, move limit
  case  5, xmin = 0;   xmax = 10; x = 10; move = 0.02; low = -inf; upp = inf; sinit = 0.2;  sincr = 1.2; sdecr = 0.7; mu = 0;  % SLP, small move limit
  case  6, xmin = 0;   xmax = 10; x =  3; move = inf;  low = 0;    upp = inf; sinit = 0.2;  sincr = 1.2; sdecr = 0.7; mu = 0;  % CONLIN
  case  7, xmin = 0;   xmax = 10; x = 10; move = inf;  low = 0;    upp = inf; sinit = 0.2;  sincr = 1.2; sdecr = 0.7; mu = 0;  % CONLIN, different starting point
  case  8, xmin = 0;   xmax = 10; x =  3; move = inf;  low = [];   upp = [];  sinit = 0.2;  sincr = 1.2; sdecr = 0.7; mu = 0;  % MMA
  case  9, xmin = 0;   xmax = 10; x =  3; move = inf;  low = [];   upp = [];  sinit = 0.02; sincr = 1.2; sdecr = 0.7; mu = 0;  % MMA, tighter asymptotes
  case 10, xmin = 0;   xmax = 10; x =  3; move = inf;  low = -20;  upp = 30;  sinit = 0.2;  sincr = 1.2; sdecr = 0.7; mu = 0;  % MMA, fixed asymptotes -> SLP
  case 11, xmin = 0;   xmax = 10; x =  3; move = inf;  low = 0;    upp = 30;  sinit = 0.2;  sincr = 1.2; sdecr = 0.7; mu = 0;  % MMA, fixed asymptotes -> CONLIN
end

% DEFINE PLOT COLORS
darkblue = [0.000 0.447 0.741];
darkred = [0.850 0.325 0.098];
black = [0.000 0.000 0.000];
darkgray = [0.400 0.400 0.400];
lightgray = [0.750 0.750 0.750];

% PREPARE FIGURE
X = [-2:0.0001:12];
F0 = 1/3*X+1/10*(X-5).^3+5;
F = 50./X.^2-2;
fprintf('\nBlue line: objective function f0(x)\n');
fprintf('Red line: constraint f1(x)\n');
fprintf('Gray lines: box constraints xmin and xmax\n\n');
figure;
plot([-1e6 1e6],[0 0],'Color',black);
hold('on');
plot([xmin xmin],[-1e6 1e6],'LineWidth',2,'Color',lightgray);
plot([xmax xmax],[-1e6 1e6],'LineWidth',2,'Color',lightgray);
hf0 = plot(nan,nan,'Color',darkblue);
hf = plot(nan,nan,'Color',darkred);
halpha = plot(nan,nan,'Color',lightgray);
hbeta = plot(nan,nan,'Color',lightgray);
plot(X,F0,'LineWidth',2,'Color',darkblue);
plot(X,F,'LineWidth',2,'Color',darkred);
hf0 = plot(nan,nan,'Color',darkblue);
hf = plot(nan,nan,'Color',darkred);
hf0x = plot(nan,nan,'o','MarkerEdgeColor',darkblue,'MarkerFaceColor',darkblue);
hfx = plot(nan,nan,'o','MarkerEdgeColor',darkred,'MarkerFaceColor',darkred);
hx = plot(nan,nan,'o','MarkerEdgeColor',darkgray,'MarkerFaceColor',darkgray);
hxnew = plot(nan,nan,'o','MarkerEdgeColor',black,'MarkerFaceColor',black);
xlim([-2 12]);
ylim([-10 10]);
xlabel('x');
ylabel('f_i');
pause;

% OPTIMIZATION LOOP
iter = 1;
mmaparams = [];
tol = 1e-4;
change = inf;
history = struct;
while change>tol

  % OBJECTIVE FUNCTION AND CONSTRAINTS
  f0 = 1/3*x+1/10*(x-5)^3+5;
  f = 50/x^2-2;

  % SENSITIVITIES
  df0dx = 1/3+3/10*(x-5)^2;
  dfdx = -100/x^3;

  % SUPPRESS WARNING FOR NEARLY SINGULAR MATRIX INVERSION IN CASE 2
  if preset==2
    w = warning('off','MATLAB:nearlySingularMatrix');
  end

  % MMA UPDATE
  [xnew,y,z,lmult,mmaparams,subp,change] = mma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams,move,low,upp,sinit,sincr,sdecr,mu,'silent');

  % RESTORE WARNING STATE
  if preset==2
    warning(w);
  end

  % COLLECT CONVERGENCE HISTORY
  history.iter(:,iter) = iter;
  history.x(:,iter) = x;
  history.low(:,iter) = subp.low;
  history.upp(:,iter) = subp.upp;
  history.f0(:,iter) = f0;
  history.f(:,iter) = f;

  % PRINT RESULTS
  fprintf('--------------------------------------------------\n');
  fprintf('ITERATION %i\n\n',mmaparams.iter);
  fprintf('Current value x = %0.5f\n',x);
  if iter>1
    fprintf('Change = %0.5f\n\n',change);
  else
    fprintf('\n');
  end
  fprintf('Objective function f0(x) = %0.5f\n',f0);
  fprintf('Constraint f1(x) = %0.5f\n\n',f);
  fprintf('Lower asymptote L = %0.5f\n',subp.low);
  fprintf('Upper asymptote U = %0.5f\n\n',subp.upp);

  % PLOT RESULTS
  [F0,F] = mmasubp(subp,X);
  set(hf0x,'XData',x,'YData',f0);
  set(hfx,'XData',x,'YData',f);
  set(hf0,'XData',X,'YData',F0);
  set(hf,'XData',X,'YData',F);
  set(hxnew,'XData',xnew,'YData',0);
  set(hx,'XData',x,'YData',0);
  set(halpha,'XData',[subp.alpha subp.alpha],'YData',[-1e6 1e6]);
  set(hbeta,'XData',[subp.beta subp.beta],'YData',[-1e6 1e6]);
  pause;

  % CHECK CONVERGENCE CRITERION AND MOVE ON TO NEXT ITERATION
  if change>tol
    iter = iter+1;
    x = xnew;
  end
end

% PRINT FINAL MESSAGE
fprintf('--------------------------------------------------\n');
fprintf('CONVERGENCE CRITERION SATISFIED - TERMINATING\n\n');

