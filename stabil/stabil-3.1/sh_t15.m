function [Ni,dN_dxi,dN_deta] = sh_t15(xi,eta)

%SH_T15    Shape functions for a 15 node triangular element.
%
%   [Ni,dN_dxi,dN_deta] = sh_t15(xi,eta) returns the shape functions and
%   its derivatives with respect to the natural coordinates in the point (xi,eta)*
%
%   xi      natural coordinate                  (scalar)
%   eta     natural coordinate                  (scalar)
%   Ni      shape functions in point (xi,eta)   (15 * 1)
%   dN_dxi  derivative of Ni to xi              (15 * 1)
%   dN_deta derivative of Ni to eta             (15 * 1)
%
%   The nodes have the following order:
%
%   see also SH_T3, SH_T6, SH_T15.

Ni=zeros(15,1);
Ni( 1)=((eta+xi-1)*(4*eta+4*xi-1)*(4*eta+4*xi-2)*(4*eta+4*xi-3))/6;
Ni( 2)=(xi*(4*xi-1)*(4*xi-2)*(4*xi-3))/6;
Ni( 3)=(eta*(4*eta-1)*(4*eta-2)*(4*eta-3))/6;
Ni( 4)=-(8*xi*(eta+xi-1)*(4*eta+4*xi-2)*(4*eta+4*xi-3))/3;
Ni( 5)=xi*(4*xi-1)*(4*eta+4*xi-3)*(4*eta+4*xi-4);
Ni( 6)=-(8*xi*(4*xi-1)*(4*xi-2)*(eta+xi-1))/3;
Ni( 7)=(8*eta*xi*(4*xi-1)*(4*xi-2))/3;
Ni( 8)=4*eta*xi*(4*eta-1)*(4*xi-1);
Ni( 9)=(8*eta*xi*(4*eta-1)*(4*eta - 2))/3;
Ni(10)=-(8*eta*(4*eta-1)*(4*eta-2)*(eta+xi-1))/3;
Ni(11)=4*eta*(4*eta-1)*(eta+xi-1)*(4*eta+4*xi-3);
Ni(12)=-(8*eta*(eta+xi-1)*(4*eta+4*xi-2)*(4*eta+4*xi-3))/3;
Ni(13)=32*eta*xi*(eta+xi-1)*(4*eta+4*xi-3);
Ni(14)=-32*eta*xi*(4*xi-1)*(eta+xi-1);
Ni(15)=-32*eta*xi*(4*eta-1)*(eta+xi-1);

if nargout > 1
  [dN_dxi,dN_deta]=deal(zeros(15,1)); 
  dN_dxi( 1)=(128*eta^3)/3+128*eta^2*xi-80*eta^2+128*eta*xi^2-160*eta*xi+(140*eta)/3+(128*xi^3)/3-80*xi^2+(140*xi)/3-25/3;
  dN_dxi( 2)=(128*xi^3)/3-48*xi^2+(44*xi)/3-1;
  dN_dxi( 3)=0;
  dN_dxi( 4)=-(128*eta^3)/3-256*eta^2*xi+96*eta^2-384*eta*xi^2+384*eta*xi-(208*eta)/3-(512*xi^3)/3+288*xi^2-(416*xi)/3+16;
  dN_dxi( 5)=128*eta^2*xi-16*eta^2+384*eta*xi^2-288*eta*xi+28*eta+256*xi^3-384*xi^2+152*xi-12;
  dN_dxi( 6)=64*eta*xi-(224*xi)/3-(16*eta)/3-128*eta*xi^2+224*xi^2-(512*xi^3)/3+16/3;
  dN_dxi( 7)=(16*eta*(24*xi^2-12*xi+1))/3;
  dN_dxi( 8)=4*eta*(4*eta-1)*(8*xi-1);
  dN_dxi( 9)=(8*eta*(4*eta-1)*(4*eta-2))/3;
  dN_dxi(10)=-(8*eta*(4*eta-1)*(4*eta-2))/3;
  dN_dxi(11)=4*eta*(4*eta-1)*(8*eta+8*xi-7);
  dN_dxi(12)=-(16*eta*(24*eta^2+48*eta*xi-36*eta+24*xi^2-36*xi+13))/3;
  dN_dxi(13)=32*eta*(4*eta^2+16*eta*xi-7*eta+12*xi^2-14*xi+3);
  dN_dxi(14)=-32*eta*(8*eta*xi-10*xi-eta+12*xi^2+1);
  dN_dxi(15)=-32*eta*(4*eta-1)*(eta+2*xi-1);

  dN_deta( 1)=(128*eta^3)/3+128*eta^2*xi-80*eta^2+128*eta*xi^2-160*eta*xi+(140*eta)/3+(128*xi^3)/3-80*xi^2+(140*xi)/3-25/3;
  dN_deta( 2)=0;
  dN_deta( 3)=(128*eta^3)/3-48*eta^2+(44*eta)/3-1;
  dN_deta( 4)=-(16*xi*(24*eta^2+48*eta*xi-36*eta+24*xi^2-36*xi+13))/3;
  dN_deta( 5)=4*xi*(4*xi-1)*(8*eta+8*xi-7);
  dN_deta( 6)=-(8*xi*(4*xi-1)*(4*xi-2))/3;
  dN_deta( 7)=(8*xi*(4*xi-1)*(4*xi-2))/3;
  dN_deta( 8)=4*xi*(8*eta-1)*(4*xi-1);
  dN_deta( 9)=(16*xi*(24*eta^2-12*eta+1))/3;
  dN_deta(10)=64*eta*xi-(16*xi)/3-(224*eta)/3-128*eta^2*xi+224*eta^2-(512*eta^3)/3+16/3;
  dN_deta(11)=256*eta^3+384*eta^2*xi-384*eta^2+128*eta*xi^2-288*eta*xi+152*eta-16*xi^2+28*xi-12;
  dN_deta(12)=-(512*eta^3)/3-384*eta^2*xi+288*eta^2-256*eta*xi^2+384*eta*xi-(416*eta)/3-(128*xi^3)/3+96*xi^2-(208*xi)/3+16;
  dN_deta(13)=32*xi*(12*eta^2+16*eta*xi-14*eta+4*xi^2-7*xi+3);
  dN_deta(14)=-32*xi*(4*xi-1)*(2*eta+xi-1);
  dN_deta(15)=-32*xi*(8*eta*xi-xi-10*eta+12*eta^2+1);
end