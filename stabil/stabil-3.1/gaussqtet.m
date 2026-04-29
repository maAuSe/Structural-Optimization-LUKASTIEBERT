function [x,H] = gaussqtet(n)

%GAUSSQTET   Gauss points for 3D numerical integration on a tetrahedron.
%
%   [x,H] = gaussqtet(n) returns the coordinates and weights for a 3D
%   gauss-legendre quadrature on a tetrahedron. The coordinates are
%   returned in natural coordinates (i.e. no volume coordinates) 
%
%   n   number of integration points 
%   x   coordinates of the integration points  (n * 3)
%   H   weights used in summation              (1 * n)
%
%   See also GAUSSQ, GAUSSQTRI 

if (n==1)
  x=[0.25 0.25 0.25];
  H=[0.166666666666667];
elseif (n==4)
  x=[0.138196601125011 0.138196601125011 0.138196601125011
     0.585410196624969 0.138196601125011 0.138196601125011
     0.138196601125011 0.585410196624969 0.138196601125011
     0.138196601125011 0.138196601125011 0.585410196624969];
  H=[0.041666666666667
     0.041666666666667
     0.041666666666667
     0.041666666666667];
elseif (n==5)
  x=[0.250000000000000 0.250000000000000 0.250000000000000
     0.500000000000000 0.166666666666667 0.166666666666667
     0.166666666666667 0.500000000000000 0.166666666666667 
     0.166666666666667 0.166666666666667 0.500000000000000
     0.166666666666667 0.166666666666667 0.166666666666667];
  H=[-0.133333333333333
      0.075000000000000
      0.075000000000000
      0.075000000000000
      0.075000000000000];     
else
  error(['Integration scheme with ' num2str(n) ' points for tetrahedral elements not available.']);  
end