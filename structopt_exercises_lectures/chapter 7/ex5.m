%EX5
%
%   B-spline surface [x(u,v),y(u,v),z(u,v)] with 3*3 control vertices.

% B-SPLINE BASIS FUNCTIONS
pu = 2;                          % Spline degree in terms of u
pv = 2;                          % spline degree in terms of v
Nu = 3;                          % Number of control vertices in u-direction
Nv = 3;                          % Number of control vertices in v-direction
u = linspace(0,1,25);            % Free parameter sampling in u-direction
v = linspace(0,1,25);            % Free parameter sampling in v-direction
Bu = bspline(pu,Nu,u);           % B-spline basis functions in u-direction
Bv = bspline(pv,Nv,v);           % B-spline basis functions in v-direction

% CONTROL VERTICES
X = [ 9.99   0.00  -9.99;
     10.61   0.00 -10.61;
      9.99   0.00  -9.99];
Y = [-2.12  -2.12  -2.12;
      0.00   0.00   0.00;
      2.12   2.12   2.12];
Z = [ 9.99  14.12   9.99;
     10.61  15.00  10.61;
      9.99  14.12   9.99];

% GENERATE SPLINE SURFACE
x = Bu*X*Bv';
y = Bu*Y*Bv';
z = Bu*Z*Bv';

% PLOT SPLINE SURFACE
figure;
mesh(x,y,z,'EdgeColor',[0 0.447 0.741]);
hold('on');
plot3(X(:),Y(:),Z(:),'o','Color',[0.850 0.325 0.098]);
axis('equal');
axis([-12 12 -6 6 9 20]);
xlabel('x');
ylabel('y');
zlabel('z');
