function [Ni,dN_dxi,dN_deta] = sh_t6(xi,eta)

%SH_T6     Shape functions for an 6 node triangular element.
%
%   [Ni,dN_dxi,dN_deta] = sh_t6(xi,eta) returns the shape functions and
%   its derivatives with respect to the natural coordinates in the point (xi,eta).
%
%   xi      natural coordinate                  (scalar)
%   eta     natural coordinate                  (scalar)
%   Ni      shape functions in point (xi,eta)   (6 * 1)
%   dN_dxi  derivative of Ni to xi              (6 * 1)
%   dN_deta derivative of Ni to eta             (6 * 1)
%
%   see also KE_SHELL6.

Ni=zeros(6,1);
Ni(1)=(1-xi-eta).*(1-2*xi-2*eta);
Ni(2)=-xi.*(1-2*xi);
Ni(3)=-eta.*(1-2*eta);
Ni(4)=4.0*xi.*(1-xi-eta);
Ni(5)=4.0*xi.*eta;
Ni(6)=4.0*eta.*(1-xi-eta);

if nargout > 1
  [dN_dxi,dN_deta]=deal(zeros(6,1));
  dN_dxi(1)=-3+4*xi+4*eta;
  dN_dxi(2)=-1+4*xi;
  dN_dxi(3)= 0.0;
  dN_dxi(4)= 4-8*xi-4*eta;
  dN_dxi(5)= 4*eta;
  dN_dxi(6)=-4*eta;
  
  dN_deta(1)=-3+4*xi+4*eta;
  dN_deta(2)= 0.0;
  dN_deta(3)=-1+4*eta;
  dN_deta(4)=-4*xi;
  dN_deta(5)= 4*xi;
  dN_deta(6)= 4-4*xi-8*eta;
end