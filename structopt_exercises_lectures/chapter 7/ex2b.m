%EX2B
%
%   B-spline basis functions with uniform knot vectors.

% FREE PARAMETER SAMPLING
u = linspace(-0.2,1.2,10000);

% PLOT 3-RD DEGREE B-SPLINE BASIS FUNCTIONS
p = 3;
for N = [4,5,6,10,20]
  B = bspline(p,N,u);
  figure;
  plot(u,B);
  ylim([-0.2 1.2]);
  xlabel('u');
  ylabel(sprintf('B_{i,%i}',p));
  title(sprintf('%i B-spline basis functions of degree %i',N,p));
end
