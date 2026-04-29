%EX2A
%
%   Exercise 4.6: use of SLP, CONLIN, and MMA to solve a four-bar truss
%   optimization problem.

% Mattias Schevenels
% February 2021

% PROBLEM FORMULATION
xmin = 0.2;       % Lower bound for x
xmax = 2.5;       % Upper bound for x
x = [2;1.5];      % Initial value

% MMA PARAMETERS
move = inf;       % Disable move limit
low = 0;          % Lower asymptote: use -INF for SLP,   0 for CONLIN, [] for MMA
upp = inf;        % Upper asymptote: use  INF for SLP, INF for CONLIN, [] for MMA
sinit = 0.5;      % Asymptotes relative initial position.
sincr = 1.2;      % Factor used to increase the distance between the asymptotes
sdecr = 0.7;      % Factor used to decrease the distance between the asymptotes
mu = 0;           % Relative distance to keep away from the asymptotes
tol = 1e-4;       % Convergence threshold

% DEFINE PLOT COLORS
darkblue = [0.000 0.447 0.741];
darkred = [0.850 0.325 0.098];
black = [0.000 0.000 0.000];
darkgray = [0.400 0.400 0.400];
lightgray = [0.750 0.750 0.750];

% PLOT EXACT OBJECTIVE FUNCTION AND CONSTRAINT
[X1,X2] = ndgrid(0:0.01:3,0:0.01:3);
F0 = X1+X2;
F = 8./(16*X1+9*X2)-4.5./(9*X1+16*X2)-0.1;
figure;
h3 = plot([xmin xmax xmax xmin xmin],[xmin xmin xmax xmax xmin],'LineWidth',2,'Color',lightgray);
hold('on');
axis('equal');
xlim([0 3]);
ylim([0 3]);
[~,h1] = contour(X1,X2,F0,[-10:1:10],'LineWidth',1.5,'Color',darkblue);
[~,h2] = contour(X1,X2,F,[0 0],'LineWidth',2,'Color',darkred);
htop = [h1;h2;h3];
xlabel('x_1');
ylabel('x_2');
fprintf('\nBlue lines: objective function f0(x)\n');
fprintf('Red line: constraint f1(x)\n');
fprintf('Gray lines: box constraints xmin and xmax\n\n');
pause;

% OPTIMIZATION LOOP
iter = 1;
change = inf;
mmaparams = [];
while change>tol

  % OBJECTIVE FUNCTION AND CONSTRAINTS
  f0 = sum(x);
  f = 8/(16*x(1)+9*x(2))-4.5/(9*x(1)+16*x(2))-0.1;

  % OBJECTIVE FUNCTION SENSITIVITIES
  df0dx1 = 1;
  df0dx2 = 1;
  df0dx  = [df0dx1;
            df0dx2];

  % CONSTRAINT SENSITIVITIES
  df1dx1 = -8/(16*x(1)+9*x(2))^2*16 + 4.5/(9*x(1)+16*x(2))^2*9;
  df1dx2 = -8/(16*x(1)+9*x(2))^2*9  + 4.5/(9*x(1)+16*x(2))^2*16;
  dfdx = [df1dx1 df1dx2];

  % MMA UPDATE
  [xnew,y,z,lmult,mmaparams,subp,change] = mma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams,move,low,upp,sinit,sincr,sdecr,mu,'silent');

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
  fprintf('Current value x = [%0.5f %0.5f]\n',x);
  if iter>1
    fprintf('Change = %0.5f\n\n',change);
  else
    fprintf('\n');
  end
  fprintf('Objective function f0(x) = %0.5f\n',f0);
  fprintf('Constraint f1(x) = %0.5f\n\n',f);
  fprintf('Lower asymptotes L = [%0.5f %0.5f]\n',subp.low);
  fprintf('Upper asymptotes U = [%0.5f %0.5f]\n\n',subp.upp);

  % PLOT RESULTS
  if iter>1
    delete(hf0);
    delete(hxnew);
    delete(hf);
    delete(halphabeta);
    delete(hx);
  end
  [F0,F] = mmasubp(subp,X1,X2);
  halphabeta = plot([subp.alpha(1) subp.beta(1) subp.beta(1) subp.alpha(1) subp.alpha(1)], ...
    [subp.alpha(2) subp.alpha(2) subp.beta(2) subp.beta(2) subp.alpha(2)],'Color',lightgray);
  [~,hf0] = contour(X1,X2,F0,[-10:2:18],'LineWidth',0.1,'Color',darkblue);
  [~,hf] = contour(X1,X2,F,[0 0],'Color',darkred);
  set(gca,'Children',[htop;setdiff(get(gca,'Children'),htop,'stable')]);
  hx = plot(x(1),x(2),'o','MarkerEdgeColor',darkgray,'MarkerFaceColor',darkgray);
  hxnew = plot(xnew(1),xnew(2),'o','MarkerEdgeColor',black,'MarkerFaceColor',black);
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

