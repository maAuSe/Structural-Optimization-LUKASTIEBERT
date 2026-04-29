function Ni = sh_t(L,b,c)

%SH_T     Shape functions for a triangular plate element.
%
%   [Ni] = sh_t(L,b,c) returns the shape functions and
%   its derivatives in point L.
%
%   L       Area coordinates [L1,L2,L3]                       (3 * 1)
%   b       Geometrical property of the triangle (see ke_dkt) (3 * 1)
%   c       Geometrical property of the triangle              (3 * 1)
%   Ni      Shape functions and derivatives                   (9 * 1)
%           in point L   
%
%   These shape functions are used to determine the mass matrix of a
%   triangular plate element.
%
%   See also KE_DKT.

Ni = zeros(1,9);
% Ni(1:3:7) = 3*L.^2-2*L.^3;
% Ni(2:3:8) = L.^2.*(b([2 3 1]).*L([3 1 2])-b([3 1 2]).*L([2 3 1]))+(b([2 3 1])-b([3 1 2]))*L(1)*L(2)*L(3);
% Ni(3:3:9) = L.^2.*(c([2 3 1]).*L([3 1 2])-c([3 1 2]).*L([2 3 1]))+(c([2 3 1])-c([3 1 2]))*L(1)*L(2)*L(3);
Ni(1:3:7) = L+L.^2.*L([2 3 1])+L.^2.*L([3 1 2])-L.*L([2 3 1]).^2-L.*L([3 1 2]).^2;
Ni(2:3:8) = b([2 3 1]).*(L([3 1 2]).*L.^2+L(1)*L(2)*L(3)/2)-b([3 1 2]).*(L.^2.*L([2 3 1])+L(1)*L(2)*L(3)/2);
Ni(3:3:9) = c([2 3 1]).*(L([3 1 2]).*L.^2+L(1)*L(2)*L(3)/2)-c([3 1 2]).*(L.^2.*L([2 3 1])+L(1)*L(2)*L(3)/2);