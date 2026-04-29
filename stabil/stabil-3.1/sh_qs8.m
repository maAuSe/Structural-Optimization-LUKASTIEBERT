function [Ni,dN_dxi,dN_deta] = sh_qs8(xi,eta)

%SH_QS8     Shape functions for an 8 node quadrilateral serendipity element.
%
%   [Ni,dN_dxi,dN_deta] = sh_qs8(xi,eta) returns the shape functions and
%   its derivatives with respect to the natural coordinates in the point (xi,eta).
%
%   xi      natural coordinate                  (scalar)
%   eta     natural coordinate                  (scalar)
%   Ni      shape functions in point (xi,eta)   (8 * 1)
%   dN_dxi  derivative of Ni to xi              (8 * 1)
%   dN_deta derivative of Ni to eta             (8 * 1)
%
%   see also KE_SHELL8.

Nodes = [ -1 -1;1  -1;1 1; -1 1;0 -1; 1 0; 0 1; -1 0];

Ni=zeros(8,1);
Ni(1:4) = 1/4*(1 + Nodes(1:4,1)*xi).*(1 + Nodes(1:4,2)*eta).*(Nodes(1:4,1)*xi+Nodes(1:4,2)*eta-1);
Ni(5) = 1/2*(1 - xi^2)*(1 + Nodes(5,2)*eta);
Ni(7) = 1/2*(1 - xi^2)*(1 + Nodes(7,2)*eta);
Ni(6) = 1/2*(1 + Nodes(6,1)*xi)*(1 - eta^2);
Ni(8) = 1/2*(1 + Nodes(8,1)*xi)*(1 - eta^2);

if nargout > 1
  [dN_dxi,dN_deta]=deal(zeros(8,1)); 
  dN_dxi(1:4,1) = 1/4*sign(Nodes(1:4,1)).*(Nodes(1:4,2)*eta + 2*Nodes(1:4,1)*xi).*(1 + Nodes(1:4,2)*eta);
  dN_dxi(5,1) = -xi*(1 + Nodes(5,2)*eta);
  dN_dxi(7,1) = -xi*(1 + Nodes(7,2)*eta);
  dN_dxi(6,1) = 1/2*sign(Nodes(6,1))*(1 - eta^2);
  dN_dxi(8,1) = 1/2*sign(Nodes(8,1))*(1 - eta^2);

  dN_deta(1:4,1) = 1/4*sign(Nodes(1:4,2)).*(2*Nodes(1:4,2)*eta + Nodes(1:4,1)*xi).*(1 + Nodes(1:4,1)*xi);
  dN_deta(5,1) = 1/2*sign(Nodes(5,2))*(1 - xi^2);
  dN_deta(7,1) = 1/2*sign(Nodes(7,2))*(1 - xi^2);
  dN_deta(6,1) = -eta*(1 + Nodes(6,1)*xi);
  dN_deta(8,1) = -eta*(1 + Nodes(8,1)*xi);
end
