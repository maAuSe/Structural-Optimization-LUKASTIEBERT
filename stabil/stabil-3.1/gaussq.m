function [x,H] = gaussq(n)

%GAUSSQ   Gauss points for 2D numerical integration.
%
%   [x,H] = gaussq(n) returns the coordinates and weights for a 2D
%   gauss-legendre quadrature.
%
%   n   number of points in one direction = 2 or 3
%   x   coordinates of gauss-points         (n^2 * 2)
%   H   weights used in summation           (1 * n^2)
%
%   See also KE_SHELL8, KELCS_SHELL4

if n == 2
    x = [-0.577350269189626 -0.577350269189626;
         0.577350269189626 -0.577350269189626;
          0.577350269189626  0.577350269189626;
         -0.577350269189626 0.577350269189626];
 
    H = [1 1 1 1]; 
end

if n == 3
    a = 0.774596669241483;
    h1 = 0.555555555555556;
    h2 = 0.888888888888889;
    x = [0 0;
        a 0;
        a a;
        0 a;
        -a 0;
        -a -a;
        0 -a;
        -a a;
        a -a];
    
    H = [h2^2 h1*h2 h1^2 h1*h2 h1*h2 h1^2 h1*h2 h1^2 h1^2];
end
end