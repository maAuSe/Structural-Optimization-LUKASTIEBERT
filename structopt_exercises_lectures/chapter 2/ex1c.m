%EX1C
%
%   Section 2.1: weight minimization of a two-bar truss subject to stress
%   constraints - solution with particle swarm optimization.

% Mattias Schevenels
% February 2017

% PROBLEM FORMULATION
n = 2;                                                      % Number of design variables
m = 2;                                                      % Number of constraints
p = 100;                                                    % Number of particles in the swarm
xmin = 1e-4;                                                % Minimum value design variables
xmax = 0.1;                                                 % Maximum value design variables
x = repmat(xmin,1,p) + rand(n,p).*repmat(xmax-xmin,1,p);    % Initial designs
rho = 7850;                                                 % Density
L = 1;                                                      % Bar length
sigma0 = 235e6;                                             % Maximum allowable stress
F = 200e3;                                                  % Force magnitude
alpha = 30/180*pi;                                          % Force angle

% OPTIMIZATION LOOP
psoparams = [];                                             % Internal PSO parameters
for iter = 1:50

  % OBJECTIVE FUNCTION AND CONSTRAINTS
  f0 = zeros(1,p);                                          % Objective function for all particles
  f = zeros(m,p);                                           % Constraints for all particles
  for k = 1:p
    f0(k) = rho*L*sum(x(:,k));                              % Total weight
    f1 = F*cos(alpha)/x(1,k)/sigma0-1;                      % Stress constraint for bar 1
    f2 = F*sin(alpha)/x(2,k)/sigma0-1;                      % Stress constraint for bar 2
    f(:,k) = [f1; f2];
  end

  % PSO UPDATE
  [x,xgbest,f0gbest,fgbest,psoparams,change,history] = pso(x,xmin,xmax,f0,f,psoparams);

end
