%EX1B
%
%   Section 2.1: weight minimization of a two-bar truss subject to stress
%   constraints - solution with MMA.

% Mattias Schevenels
% February 2017

% PROBLEM FORMULATION
xmin = 1e-4;                                      % Minimum value design variables
xmax = 0.1;                                       % Maximum value design variables
x = [0.001; 0.001];                               % Initial design
rho = 7850;                                       % Density
L = 1;                                            % Bar length
sigma0 = 235e6;                                   % Maximum allowable stress
F = 200e3;                                        % Force magnitude
alpha = 30/180*pi;                                % Force angle

% OPTIMIZATION LOOP
iter = 1;                                         % Iteration counter
change = inf;                                     % (Fictitious) initial design change
tol = 1e-7;                                       % Convergence threshold
mmaparams = [];                                   % Internal MMA parameters
while change>tol

  % OBJECTIVE FUNCTION
  f0 = rho*L*sum(x);                              % Total weight

  % CONSTRAINTS
  f1 = F*cos(alpha)/x(1)/sigma0-1;                % Stress constraint for bar 1
  f2 = F*sin(alpha)/x(2)/sigma0-1;                % Stress constraint for bar 2
  f = [f1; f2];

  % OBJECTIVE FUNCTION SENSITIVITIES
  df0dx1 = rho*L;
  df0dx2 = rho*L;
  df0dx  = [df0dx1;
            df0dx2];

  % CONSTRAINT SENSITIVITIES
  df1dx1 = -F*cos(alpha)/x(1)^2/sigma0;
  df1dx2 = 0;
  df2dx1 = 0;
  df2dx2 = -F*sin(alpha)/x(2)^2/sigma0;
  dfdx = [df1dx1 df1dx2;
          df2dx1 df2dx2];

  % MMA UPDATE
  [xnew,y,z,lmult,mmaparams,subp,change,history] = mma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams);

  % CHECK CONVERGENCE CRITERION AND MOVE ON TO NEXT ITERATION
  if change>tol
    iter = iter+1;
    x = xnew;
  end

end
