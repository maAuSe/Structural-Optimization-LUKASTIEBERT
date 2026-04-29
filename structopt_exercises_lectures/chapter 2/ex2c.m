%EX2C
%
%   Section 2.2: weight minimization of a two-bar truss subject to stress
%   and instability constraints - solution with particle swarm optimization.

%   This code contains the solution of an exercise introduced in the web lectures.
%   Try to solve the exercise independently before inspecting or running the code.

% Mattias Schevenels
% February 2017

% PROBLEM FORMULATION
n = 2;                                                      % Number of design variables
m = 3;                                                      % Number of constraints
p = 100;                                                    % Number of particles in the swarm
xmin = 1e-4;                                                % Minimum value design variables
xmax = 0.1;                                                 % Maximum value design variables
x = repmat(xmin,1,p) + rand(n,p).*repmat(xmax-xmin,1,p);    % Initial designs
rho = 7850;                                                 % Density
L = 1;                                                      % Bar length
E = 210e9;                                                  % Young's modulus
sigma0 = 235e6;                                             % Maximum allowable stress
F = 200e3;                                                  % Force magnitude

% OPTIMIZATION LOOP
psoparams = [];                                             % Internal PSO parameters
for iter = 1:50

  % OBJECTIVE FUNCTION AND CONSTRAINTS
  f0 = zeros(1,p);                                          % Objective function for all particles
  f = zeros(m,p);                                           % Constraints for all particles
  for k = 1:p
    f0(k) = rho*L*sum(x(:,k));                              % Total weight
    f1 = F/sqrt(2)/x(1,k)/sigma0-1;                         % Stress constraint for bar 1
    f2 = F/sqrt(2)/x(2,k)/sigma0-1;                         % Stress constraint for bar 2
    f3 = 16*F*L^2/sqrt(2)/pi/E/x(2,k)^2-1;                  % Instability constraint for bar 2
    f(:,k) = [f1; f2; f3];
  end

  % PSO UPDATE
  [x,xgbest,f0gbest,fgbest,psoparams,change,history] = pso(x,xmin,xmax,f0,f,psoparams);

end
