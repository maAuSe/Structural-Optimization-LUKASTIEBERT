%EX3C
%
%   Section 2.4: weight minimization of a cantilever beam subject to a
%   displacement constraint - solution with particle swarm optimization.

% Mattias Schevenels
% February 2017

% PROBLEM FORMULATION
n = 2;                                                              % Number of design variables (beam segments)
m = 1;                                                              % Number of constraints
p = 100;                                                            % Number of particles in the swarm
xmin = 1e-4;                                                        % Minimum value design variables
xmax = 1;                                                           % Maximum value design variables
x = repmat(xmin,1,p) + rand(n,p).*repmat(xmax-xmin,1,p);            % Initial designs
L = 1/n;                                                            % Segment length
rho = 7850;                                                         % Density
E = 210e9;                                                          % Young's modulus
t = 0.01;                                                           % Wall thickness
F = 10e3;                                                           % Force
delta0 = 0.004;                                                     % Maximum allowable displacement

% OPTIMIZATION LOOP
psoparams = [];                                                     % Internal PSO parameters
for iter = 1:50

  % OBJECTIVE FUNCTION AND CONSTRAINTS
  f0 = zeros(1,p);                                                  % Objective function for all particles
  f = zeros(m,p);                                                   % Constraints for all particles
  for k = 1:p
    f0(k) = 4*rho*L*t*sum(x(:,k));                                  % Total weight
    A = [1:n]';                                                     % Segment counter
    f(:,k) = 3*F*L^3/2/E/t*sum((A.^2-A+1/3).*x(:,k).^-3)/delta0-1;  % Beam tip displacement
  end

  % PSO UPDATE
  [x,xgbest,f0gbest,fgbest,psoparams,change,history] = pso(x,xmin,xmax,f0,f,psoparams);

  % PLOT CURRENT DESIGN
  plot(reshape([1:-L:L;1-L:-L:0],[],1),reshape([xgbest,xgbest]',[],1));
  ylim([0 0.2]);
  drawnow;

end
