%EX3B
%
%   Section 2.4: weight minimization of a cantilever beam subject to a
%   displacement constraint - solution with MMA.

% Mattias Schevenels
% February 2017

% PROBLEM FORMULATION
n = 2;                                                  % Number of beam segments
xmin = 1e-4;                                            % Minimum value design variables
xmax = 1;                                               % Maximum value design variables
x = repmat(0.1,n,1);                                    % Initial design
L = 1/n;                                                % Segment length
rho = 7850;                                             % Density
E = 210e9;                                              % Young's modulus
t = 0.01;                                               % Wall thickness
F = 10e3;                                               % Force
delta0 = 0.004;                                         % Maximum allowable displacement

% OPTIMIZATION LOOP
iter = 1;                                               % Iteration counter
change = inf;                                           % Change in design variables wrt previous iteration
tol = 1e-4;                                             % Convergence threshold
mmaparams = [];                                         % Internal MMA parameters
while change>tol

  % OBJECTIVE FUNCTION
  f0 = 4*rho*L*t*sum(x);                                % Total weight

  % CONSTRAINTS
  A = [1:n]';                                           % Segment counter
  f = 3*F*L^3/2/E/t*sum((A.^2-A+1/3).*x.^-3)/delta0-1;  % Beam tip displacement

  % OBJECTIVE FUNCTION SENSITIVITIES
  df0dx = repmat(4*rho*L*t,n,1);

  % CONSTRAINT SENSITIVITIES
  dfdx = -9*F*L^3/2/E/t*((A.^2-A+1/3).*x.^-4)/delta0;

  % MMA UPDATE
  [xnew,y,z,lmult,mmaparams,subp,change,history] = mma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams);

  % PLOT CURRENT DESIGN
  plot(reshape([1:-L:L;1-L:-L:0],[],1),reshape([x,x]',[],1));
  ylim([0 0.2]);
  drawnow;

  % CHECK CONVERGENCE CRITERION AND MOVE ON TO NEXT ITERATION
  if change>tol
    iter = iter+1;
    x = xnew;
  end

end
