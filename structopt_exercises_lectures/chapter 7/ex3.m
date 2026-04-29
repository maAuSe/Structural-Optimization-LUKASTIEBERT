%EX3
%
%   B-spline x(u) with 5 control vertices.

% B-SPLINE BASIS FUNCTIONS
p = 2;                           % Spline degree
N = 5;                           % Number of control vertices
u = linspace(-0.2,1.2,1000);     % Free parameter sampling
B = bspline(p,N,u);

% CONTROL VERTICES
X = [1 3 4 3 5]';

% GENERATE SPLINE
x = B*X;

% PLOT SPLINE
figure;
plot(u,B*diag(X));
hold('on');
plot(u,x,'k','LineWidth',2);
ylim([-0.2 6]);
xlabel('u');
ylabel('x');
