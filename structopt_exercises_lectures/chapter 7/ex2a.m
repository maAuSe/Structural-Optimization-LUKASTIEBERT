%EX2A
%
%   Bezier spline basis functions.

% FREE PARAMETER SAMPLING
u = linspace(-0.2,1.2,10000);

% PLOT 0-TH TO 5-RD DEGREE BEZIER SPLINE BASIS FUNCTIONS
for p = 0:5
  B = bspline(p,[],u);
  figure;
  plot(u,B);
  ylim([-0.2 1.2]);
  xlabel('u');
  ylabel(sprintf('B_{i,%i}',p));
  title(sprintf('Bezier spline basis functions of degree %i',p));
end
