%EX1A
%
%   Section 2.1: weight minimization of a two-bar truss subject to stress
%   constraints - analytical solution.

% Mattias Schevenels
% February 2017

% PROBLEM FORMULATION
rho = 7850;                   % Density
L = 1;                        % Bar length
sigma0 = 235e6;               % Maximum allowable stress
F = 200e3;                    % Force magnitude
alpha = 30/180*pi;            % Force angle

% ANALYTICAL SOLUTION
A1 = F*cos(alpha)/sigma0;     % Bar 1 section
A2 = F*sin(alpha)/sigma0;     % Bar 2 section
W = rho*L*(A1+A2);            % Total weight

% VERIFY CONSTRAINTS
sigma1 = F/A1*cos(alpha);     % Bar 1 stress
sigma2 = F/A2*sin(alpha);     % Bar 2 stress
