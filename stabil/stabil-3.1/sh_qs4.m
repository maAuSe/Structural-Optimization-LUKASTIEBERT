function [Ni,dN_dxi,dN_deta] = sh_qs4(xi,eta)

%SH_QS4     Shape functions for a quadrilateral serendipity element with 4
%           nodes.
%
%   [Ni,dN_dxi,dN_deta] = sh_qs4(xi,eta) returns the shape functions and
%   its derivatives to the natural coordinates in the point (xi,eta).
%
%   xi      Natural coordinate                  (scalar)
%   eta     Natural coordinate                  (scalar)
%   Ni      Shape functions in point (xi,eta)   (4 * 1)
%   dN_dxi  Derivative of Ni to xi              (4 * 1)
%   dN_deta Derivative of Ni to eta             (4 * 1)
%
%   See also KELCS_SHELL4.


Nodes = [ -1 -1;1  -1;1 1; -1 1];

Ni = 1/4*(1 + Nodes(:,1)*xi).*(1 + Nodes(:,2)*eta);

if nargout > 1, dN_dxi = 1/4*Nodes(:,1).*(1 + Nodes(:,2)*eta); end
if nargout > 2, dN_deta = 1/4*Nodes(:,2).*(1 + Nodes(:,1)*xi); end

end