%EX4
%
%   B-spline [x(u),y(u)] with 8 control vertices.

% B-SPLINE BASIS FUNCTIONS
p = 3;                           % Spline degree
N = 8;                           % Number of control vertices
u = linspace(0,1,5000);          % Free parameter sampling
B = bspline(p,N,u);

% CONTROL VERTICES
X = [0 1 2 5 5 6 7 8]';
Y = [0 2 1 1 4 6 3 0]';

% GENERATE SPLINE
x = B*X;
y = B*Y;

% PLOT SPLINE
figure;
plot(x,y,X,Y,':o');
xlim([-1 9]);
ylim([-1 7]);
xlabel('x');
ylabel('y');
