function [Ni,dNi_dxi,dNi_deta,dNi_dzeta] = sh_vs20(xi,eta,zeta)

%SH_VS20    Shape functions for a volume serendipity element with 20 nodes.
%
%   [Ni,dNi_dxi,dNi_deta,dNi_dzeta] = sh_vs20(xi,eta,zeta) 
%   returns the shape functions and its derivatives to the natural coordinates
%   in the point (xi,eta,zeta).
%
%   xi       Natural coordinate                  (scalar)
%   eta      Natural coordinate                  (scalar)
%   zeta     Natural coordinate                  (scalar)
%   Ni       Shape functions in point (xi,eta)   (20 * 1)
%   dN_dxi   Derivative of Ni to xi              (20 * 1)
%   dN_deta  Derivative of Ni to eta             (20 * 1)
%   dN_dzeta Derivative of Ni to zeta            (20 * 1)
%
%   See also KE_SOLID20.

% --- Shape function ---
Ni=zeros(20,1);
Ni( 1)= 0.125.*((xi-1).*(eta-1).*(zeta-1).*(xi+eta+zeta+2));  % N01(xi1,xi2,xi3)
Ni( 2)=-0.125.*((xi+1).*(eta-1).*(zeta-1).*(eta-xi+zeta+2));  % N02(xi1,xi2,xi3)
Ni( 3)=-0.125.*((xi+1).*(eta+1).*(zeta-1).*(xi+eta-zeta-2));  % N03(xi1,xi2,xi3)
Ni( 4)=-0.125.*((xi-1).*(eta+1).*(zeta-1).*(xi-eta+zeta+2));  % N04(xi1,xi2,xi3)
Ni( 5)=-0.125.*((xi-1).*(eta-1).*(zeta+1).*(xi+eta-zeta+2));  % N05(xi1,xi2,xi3)
Ni( 6)=-0.125.*((xi+1).*(eta-1).*(zeta+1).*(xi-eta+zeta-2));  % N06(xi1,xi2,xi3)
Ni( 7)= 0.125.*((xi+1).*(eta+1).*(zeta+1).*(xi+eta+zeta-2));  % N07(xi1,xi2,xi3)
Ni( 8)= 0.125.*((xi-1).*(eta+1).*(zeta+1).*(xi-eta-zeta+2));  % N08(xi1,xi2,xi3)
Ni( 9)=-0.25.*(xi.^2 - 1).*(eta - 1)   .*(zeta - 1);          % N09(xi1,xi2,xi3)
Ni(10)= 0.25.*(xi + 1)   .*(eta.^2 - 1).*(zeta - 1);          % N10(xi1,xi2,xi3)
Ni(11)= 0.25.*(xi.^2 - 1).*(eta + 1)   .*(zeta - 1);          % N11(xi1,xi2,xi3)
Ni(12)=-0.25.*(xi - 1)   .*(eta.^2 - 1).*(zeta - 1);          % N12(xi1,xi2,xi3)
Ni(13)= 0.25.*(xi.^2 - 1).*(eta - 1)   .*(zeta + 1);          % N13(xi1,xi2,xi3)
Ni(14)=-0.25.*(xi + 1)   .*(eta.^2 - 1).*(zeta + 1);          % N14(xi1,xi2,xi3)
Ni(15)=-0.25.*(xi.^2 - 1).*(eta + 1)   .*(zeta + 1);          % N15(xi1,xi2,xi3)
Ni(16)= 0.25.*(xi - 1)   .*(eta.^2 - 1).*(zeta + 1);          % N16(xi1,xi2,xi3)
Ni(17)=-0.25.*(xi - 1)   .*(zeta.^2 - 1).*(eta - 1);          % N17(xi1,xi2,xi3)
Ni(18)= 0.25.*(xi + 1)   .*(zeta.^2 - 1).*(eta - 1);          % N18(xi1,xi2,xi3)
Ni(19)=-0.25.*(xi + 1)   .*(zeta.^2 - 1).*(eta + 1);          % N19(xi1,xi2,xi3)
Ni(20)= 0.25.*(xi - 1)   .*(zeta.^2 - 1).*(eta + 1);          % N20(xi1,xi2,xi3)


% --- Shape function derivative ---
if nargout > 1
  [dNi_dxi,dNi_deta,dNi_dzeta] = deal(zeros(20,1));
  dNi_dxi( 1,1)= ((eta - 1).*(zeta - 1).*(2.*xi + eta + zeta + 1))/8;  % dN01(xi1,xi2,xi3)/dxi1
  dNi_dxi( 2,1)=-((eta - 1).*(zeta - 1).*(eta - 2.*xi + zeta + 1))/8;  % dN02(xi1,xi2,xi3)/dxi1
  dNi_dxi( 3,1)=-((eta + 1).*(zeta - 1).*(2.*xi + eta - zeta - 1))/8;  % dN03(xi1,xi2,xi3)/dxi1
  dNi_dxi( 4,1)=-((eta + 1).*(zeta - 1).*(2.*xi - eta + zeta + 1))/8;  % dN04(xi1,xi2,xi3)/dxi1
  dNi_dxi( 5,1)=-((eta - 1).*(zeta + 1).*(2.*xi + eta - zeta + 1))/8;  % dN05(xi1,xi2,xi3)/dxi1
  dNi_dxi( 6,1)=-((eta - 1).*(zeta + 1).*(2.*xi - eta + zeta - 1))/8;  % dN06(xi1,xi2,xi3)/dxi1
  dNi_dxi( 7,1)= ((eta + 1).*(zeta + 1).*(2.*xi + eta + zeta - 1))/8;  % dN07(xi1,xi2,xi3)/dxi1
  dNi_dxi( 8,1)= ((eta + 1).*(zeta + 1).*(2.*xi - eta - zeta + 1))/8;  % dN08(xi1,xi2,xi3)/dxi1
  dNi_dxi( 9,1)= -(xi.*(eta - 1).*(zeta - 1))/2;                       % dN09(xi1,xi2,xi3)/dxi1
  dNi_dxi(10,1)= ((eta.^2 - 1).*(zeta - 1))/4;                         % dN10(xi1,xi2,xi3)/dxi1
  dNi_dxi(11,1)=  (xi.*(eta + 1).*(zeta - 1))/2;                       % dN11(xi1,xi2,xi3)/dxi1
  dNi_dxi(12,1)=-((eta.^2 - 1).*(zeta - 1))/4;                         % dN12(xi1,xi2,xi3)/dxi1
  dNi_dxi(13,1)=  (xi.*(eta - 1).*(zeta + 1))/2;                       % dN13(xi1,xi2,xi3)/dxi1
  dNi_dxi(14,1)=-((eta.^2 - 1).*(zeta + 1))/4;                         % dN14(xi1,xi2,xi3)/dxi1
  dNi_dxi(15,1)=-(xi.*(eta + 1).*(zeta + 1))/2;                        % dN15(xi1,xi2,xi3)/dxi1
  dNi_dxi(16,1)=((eta.^2 - 1).*(zeta + 1))/4;                          % dN16(xi1,xi2,xi3)/dxi1
  dNi_dxi(17,1)=-((zeta.^2 - 1).*(eta - 1))/4;                         % dN17(xi1,xi2,xi3)/dxi1
  dNi_dxi(18,1)=((zeta.^2 - 1).*(eta - 1))/4;                          % dN18(xi1,xi2,xi3)/dxi1
  dNi_dxi(19,1)=-((zeta.^2 - 1).*(eta + 1))/4;                         % dN19(xi1,xi2,xi3)/dxi1
  dNi_dxi(20,1)=((zeta.^2 - 1).*(eta + 1))/4;                          % dN20(xi1,xi2,xi3)/dxi1
  
  dNi_deta( 1,1)=((xi - 1).*(zeta - 1).*(xi + 2.*eta + zeta + 1))/8;   % dN01(xi1,xi2,xi3)/dxi2
  dNi_deta( 2,1)=-((xi + 1).*(zeta - 1).*(2.*eta - xi + zeta + 1))/8;  % dN02(xi1,xi2,xi3)/dxi2
  dNi_deta( 3,1)=-((xi + 1).*(zeta - 1).*(xi + 2.*eta - zeta - 1))/8;  % dN03(xi1,xi2,xi3)/dxi2
  dNi_deta( 4,1)=-((xi - 1).*(zeta - 1).*(xi - 2.*eta + zeta + 1))/8;  % dN04(xi1,xi2,xi3)/dxi2
  dNi_deta( 5,1)=-((xi - 1).*(zeta + 1).*(xi + 2.*eta - zeta + 1))/8;  % dN05(xi1,xi2,xi3)/dxi2
  dNi_deta( 6,1)=-((xi + 1).*(zeta + 1).*(xi - 2.*eta + zeta - 1))/8;  % dN06(xi1,xi2,xi3)/dxi2
  dNi_deta( 7,1)=((xi + 1).*(zeta + 1).*(xi + 2.*eta + zeta - 1))/8;   % dN07(xi1,xi2,xi3)/dxi2
  dNi_deta( 8,1)=((xi - 1).*(zeta + 1).*(xi - 2.*eta - zeta + 1))/8;   % dN08(xi1,xi2,xi3)/dxi2
  dNi_deta( 9,1)=-((xi.^2 - 1).*(zeta - 1))/4;                         % dN09(xi1,xi2,xi3)/dxi2
  dNi_deta(10,1)=(eta.*(xi + 1).*(zeta - 1))/2;                        % dN10(xi1,xi2,xi3)/dxi2
  dNi_deta(11,1)=((xi.^2 - 1).*(zeta - 1))/4;                          % dN11(xi1,xi2,xi3)/dxi2
  dNi_deta(12,1)=-(eta.*(xi - 1).*(zeta - 1))/2;                       % dN12(xi1,xi2,xi3)/dxi2
  dNi_deta(13,1)=((xi.^2 - 1).*(zeta + 1))/4;                          % dN13(xi1,xi2,xi3)/dxi2
  dNi_deta(14,1)=-(eta.*(xi + 1).*(zeta + 1))/2;                       % dN14(xi1,xi2,xi3)/dxi2
  dNi_deta(15,1)=-((xi.^2 - 1).*(zeta + 1))/4;                         % dN15(xi1,xi2,xi3)/dxi2
  dNi_deta(16,1)=(eta.*(xi - 1).*(zeta + 1))/2;                        % dN16(xi1,xi2,xi3)/dxi2
  dNi_deta(17,1)=-((zeta.^2 - 1).*(xi - 1))/4;                         % dN17(xi1,xi2,xi3)/dxi2
  dNi_deta(18,1)=((zeta.^2 - 1).*(xi + 1))/4;                          % dN18(xi1,xi2,xi3)/dxi2
  dNi_deta(19,1)=-((zeta.^2 - 1).*(xi + 1))/4;                         % dN19(xi1,xi2,xi3)/dxi2
  dNi_deta(20,1)=((zeta.^2 - 1).*(xi - 1))/4;                          % dN20(xi1,xi2,xi3)/dxi2
  
  dNi_dzeta( 1,1)=((xi - 1).*(eta - 1).*(xi + eta + 2.*zeta + 1))/8;   % dN01(xi1,xi2,xi3)/dxi3
  dNi_dzeta( 2,1)=-((xi + 1).*(eta - 1).*(eta - xi + 2.*zeta + 1))/8;  % dN02(xi1,xi2,xi3)/dxi3
  dNi_dzeta( 3,1)=-((xi + 1).*(eta + 1).*(xi + eta - 2.*zeta - 1))/8;  % dN03(xi1,xi2,xi3)/dxi3
  dNi_dzeta( 4,1)=-((xi - 1).*(eta + 1).*(xi - eta + 2.*zeta + 1))/8;  % dN04(xi1,xi2,xi3)/dxi3
  dNi_dzeta( 5,1)=-((xi - 1).*(eta - 1).*(xi + eta - 2.*zeta + 1))/8;  % dN05(xi1,xi2,xi3)/dxi3
  dNi_dzeta( 6,1)=-((xi + 1).*(eta - 1).*(xi - eta + 2.*zeta - 1))/8;  % dN06(xi1,xi2,xi3)/dxi3
  dNi_dzeta( 7,1)=((xi + 1).*(eta + 1).*(xi + eta + 2.*zeta - 1))/8;   % dN07(xi1,xi2,xi3)/dxi3
  dNi_dzeta( 8,1)=((xi - 1).*(eta + 1).*(xi - eta - 2.*zeta + 1))/8;   % dN08(xi1,xi2,xi3)/dxi3
  dNi_dzeta( 9,1)=-((xi.^2 - 1).*(eta - 1))/4;                         % dN09(xi1,xi2,xi3)/dxi3
  dNi_dzeta(10,1)=((eta.^2 - 1).*(xi + 1))/4;                          % dN10(xi1,xi2,xi3)/dxi3
  dNi_dzeta(11,1)=((xi.^2 - 1).*(eta + 1))/4;                          % dN11(xi1,xi2,xi3)/dxi3
  dNi_dzeta(12,1)=-((eta.^2 - 1).*(xi - 1))/4;                         % dN12(xi1,xi2,xi3)/dxi3
  dNi_dzeta(13,1)=((xi.^2 - 1).*(eta - 1))/4;                          % dN13(xi1,xi2,xi3)/dxi3
  dNi_dzeta(14,1)=-((eta.^2 - 1).*(xi + 1))/4;                         % dN14(xi1,xi2,xi3)/dxi3
  dNi_dzeta(15,1)=-((xi.^2 - 1).*(eta + 1))/4;                         % dN15(xi1,xi2,xi3)/dxi3
  dNi_dzeta(16,1)=((eta.^2 - 1).*(xi - 1))/4;                          % dN16(xi1,xi2,xi3)/dxi3
  dNi_dzeta(17,1)=-(zeta.*(xi - 1).*(eta - 1))/2;                      % dN17(xi1,xi2,xi3)/dxi3
  dNi_dzeta(18,1)=(zeta.*(xi + 1).*(eta - 1))/2;                       % dN18(xi1,xi2,xi3)/dxi3
  dNi_dzeta(19,1)=-(zeta.*(xi + 1).*(eta + 1))/2;                      % dN19(xi1,xi2,xi3)/dxi3
  dNi_dzeta(20,1)=(zeta.*(xi - 1).*(eta + 1))/2;                       % dN20(xi1,xi2,xi3)/dxi3
end

end