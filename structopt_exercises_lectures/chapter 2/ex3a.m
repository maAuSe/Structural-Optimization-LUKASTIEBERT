%EX3A
%
%   Section 2.4: weight minimization of a 2-segment cantilever beam subject to a
%   displacement constraint - analytical solution.

% Mattias Schevenels
% February 2017

% PROBLEM FORMULATION
L = 1/2;                                                    % Segment length
rho = 7850;                                                 % Density
E = 210e9;                                                  % Young's modulus
t = 0.01;                                                   % Wall thickness
F = 10e3;                                                   % Force
delta0 = 0.004;                                             % Maximum allowable displacement

% ANALYTICAL SOLUTION
C1 = 4*rho*L*t;                                             % Constant in objective function
C2 = 2*delta0*E*t/F/L^3;                                    % Normalized displacement
x1 = ((1+7^(1/4))/C2)^(1/3);                                % Section height for segment 1
x2 = 7^(1/4)*((1+7^(1/4))/C2)^(1/3);                        % Section height for segment 2
W = C1*(x1+x2);                                             % Total weight

% VERIFY CONSTRAINT
x = [x1; x2];                                               % Section heights
A = [1:2]';                                                 % Segment counter
delta = 3*F*L^3/2/E/t*sum((A.^2-A+1/3).*(1./x.^3));         % Beam tip displacement