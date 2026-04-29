function [Ni,dN_dxi,dN_deta] = sh_t3(xi,eta)

%SH_T3     Shape functions for a 3 node triangular element.
%
%   [Ni,dN_dxi,dN_deta] = sh_t3(xi,eta) returns the shape functions and
%   its derivatives with respect to the natural coordinates in the point (xi,eta).
%
%   xi      natural coordinate                  (scalar)
%   eta     natural coordinate                  (scalar)
%   Ni      shape functions in point (xi,eta)   (3 * 1)
%   dN_dxi  derivative of Ni to xi              (3 * 1)
%   dN_deta derivative of Ni to eta             (3 * 1)

Ni=zeros(3,1);
Ni(1,1)=(1-xi-eta);
Ni(2,1)=xi;
Ni(3,1)=eta;

if nargout > 1
  [dN_dxi,dN_deta]=deal(zeros(3,1));
  dN_dxi(1,1)=-1;
  dN_dxi(2,1)= 1;
  dN_dxi(3,1)= 0;
    
  dN_deta(1,1)=-1;
  dN_deta(2,1)= 0;
  dN_deta(3,1)= 1;
end