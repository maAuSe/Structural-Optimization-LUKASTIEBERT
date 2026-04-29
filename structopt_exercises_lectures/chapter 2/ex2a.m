%EX2A
%
%   Section 2.2: weight minimization of a two-bar truss subject to stress
%   and instability constraints - analytical solution.

%   This code contains the solution of an exercise introduced in the web lectures.
%   Try to solve the exercise independently before inspecting or running the code.

% Mattias Schevenels
% February 2017

% PROBLEM FORMULATION
rho = 7850;                                                 % Density
L = 1;                                                      % Bar length
E = 210e9;                                                  % Young's modulus
sigma0 = 235e6;                                             % Maximum allowable stress
F = 200e3;                                                  % Force magnitude

% ANALYTICAL SOLUTION
A1 = F/sqrt(2)/sigma0;                                      % Bar 1 section
A2 = max(F/sqrt(2)/sigma0,sqrt(16*F*L^2/sqrt(2)/pi/E));     % Bar 2 section
W = rho*L*(A1+A2);                                          % Total weight

% VERIFY CONSTRAINTS
sigma1 = F/sqrt(2)/A1;                                      % Bar 1 stress
sigma2 = -F/sqrt(2)/A2;                                     % Bar 2 stress
I2 = A2^2/4/pi;                                             % Bar 2 moment of inertia
Pc = pi^2*E*I2/L^2;                                         % Bar 2 buckling load
N2 = -F/sqrt(2);                                            % Bar 2 normal force
