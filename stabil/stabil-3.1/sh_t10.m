function [Ni,dN_dxi,dN_deta] = sh_t10(xi,eta)
%SH_T10    Shape functions for a 10 node triangular element.
%
%   [Ni,dN_dxi,dN_deta] = sh_t10(xi,eta) returns the shape functions and
%   its derivatives with respect to the natural coordinates in the point (xi,eta)
%
%   xi      natural coordinate                  (scalar)
%   eta     natural coordinate                  (scalar)
%   Ni      shape functions in point (xi,eta)   (10 * 1)
%   dN_dxi  derivative of Ni to xi              (10 * 1)
%   dN_deta derivative of Ni to eta             (10 * 1)
%
%   The nodes have the following order:
%
%   see also SH_T3, SH_T6, SH_T15.

Ni=zeros(10,1);
Ni( 1)=-(eta + xi - 1/3)*(eta + xi - 2/3)*((9*eta)/2 + (9*xi)/2 - 9/2);
Ni( 2)=(9*xi*(xi - 1/3)*(xi - 2/3))/2;
Ni( 3)=(9*eta*(eta - 1/3)*(eta - 2/3))/2;
Ni( 4)= xi*(eta + xi - 2/3)*((27*eta)/2 + (27*xi)/2 - 27/2);
Ni( 5)=-xi*(xi - 1/3)*((27*eta)/2 + (27*xi)/2 - 27/2);
Ni( 6)=(27*eta*xi*(xi - 1/3))/2;
Ni( 7)=(27*eta*xi*(eta - 1/3))/2;
Ni( 8)=-(27*eta*(eta - 1/3)*(eta + xi - 1))/2;
Ni( 9)= (27*eta*(eta + xi - 1)*(eta + xi - 2/3))/2;
Ni(10)=-eta*xi*(27*eta + 27*xi - 27);

if nargout > 1
  [dN_dxi,dN_deta]=deal(zeros(10,1));
  dN_dxi( 1)=18*eta + 18*xi - 27*eta*xi - (27*eta^2)/2 - (27*xi^2)/2 - 11/2;
  dN_dxi( 2)=(27*xi^2)/2 - 9*xi + 1;
  dN_dxi( 3)=0;
  dN_dxi( 4)=(27*eta^2)/2 + 54*eta*xi - (45*eta)/2 + (81*xi^2)/2 - 45*xi + 9;
  dN_dxi( 5)=(9*eta)/2 + 36*xi - 27*eta*xi - (81*xi^2)/2 - 9/2;
  dN_dxi( 6)=(9*eta*(6*xi - 1))/2;
  dN_dxi( 7)=(27*eta*(eta - 1/3))/2;
  dN_dxi( 8)=-(27*eta*(eta - 1/3))/2;
  dN_dxi( 9)=(9*eta*(6*eta + 6*xi - 5))/2;
  dN_dxi(10)=-27*eta*(eta + 2*xi - 1);

  dN_deta( 1)=18*eta + 18*xi - 27*eta*xi - (27*eta^2)/2 - (27*xi^2)/2 - 11/2;
  dN_deta( 2)=0;
  dN_deta( 3)=(27*eta^2)/2 - 9*eta + 1;
  dN_deta( 4)=(9*xi*(6*eta + 6*xi - 5))/2;
  dN_deta( 5)=-(27*xi*(xi - 1/3))/2;
  dN_deta( 6)= (27*xi*(xi - 1/3))/2;
  dN_deta( 7)=(9*xi*(6*eta - 1))/2;
  dN_deta( 8)=36*eta + (9*xi)/2 - 27*eta*xi - (81*eta^2)/2 - 9/2;
  dN_deta( 9)=(81*eta^2)/2 + 54*eta*xi - 45*eta + (27*xi^2)/2 - (45*xi)/2 + 9;
  dN_deta(10)=-27*xi*(2*eta + xi - 1);
end