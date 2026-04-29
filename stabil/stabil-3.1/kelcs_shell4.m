function [KeLCS,MeLCS] = kelcs_shell4(Node_lc,h,E,nu,rho,Options)

%KELCS_SHELL4   shell element stiffness and mass matrix in element
%               coordinate system.
%
%   [Ke,Me] = kelcs_shell4(Node_lc,h,E,nu,rho)
%    Ke     = kelcs_shell4(Node_lc,h,E,nu) 
%   returns the element stiffness and mass matrix in the element
%   coordinate system for a four node shell element (isotropic material).
%
%   Node       Node definitions           [x y z] (4 * 3)
%              Nodes should have the following order:
%              4---------3
%              |         |
%              |         |
%              |         |
%              1---------2
%
%   h          Shell thickness         
%   E          Young's modulus
%   nu         Poisson coefficient
%   rho        Mass density
%   Options    Element options [NOT available yet]    {Option1 Option2 ...}
%   KeLCS      Element stiffness matrix (24 * 24)
%   MeLCS      Element mass matrix (24 * 24)
%
%   This element is a flat shell element that consists of a bilinear
%   membrane element and four overlaid DKT triangles for the bending 
%   stiffness.
%
%   See also KE_SHELL8, ASMKM, KE_DKT.

% Miche Jansen
% 2009

K1 = zeros(8,8);
K21 = zeros(4,8);
K2 = zeros(4,4);
K12 = zeros(8,4);

KeLCS = zeros(24,24);
if nargout > 1
MeLCS = zeros(24,24);
end
%membrane stiffness
[xm,Hm]= gaussq(2);

Dm = E*h/(1-nu^2)*[1 nu 0;
     nu 1 0;
     0 0 (1-nu)/2];

J0 = [0.25*[-1 1 1 -1]*Node_lc(:,[1,2]);0.25*[-1 -1 1 1]*Node_lc(:,[1,2])];
 
for iGauss = 1:size(xm,1)
    
    xi= xm(iGauss,1);
    eta = xm(iGauss,2);
    [Ni,dN_dxi,dN_deta] = sh_qs4(xi,eta);
    
    %xg= Ni.'*Node_lc(:,1);
    %yg= Ni.'*Node_lc(:,2);
    
    J = [dN_dxi.'*Node_lc(:,1) dN_dxi.'*Node_lc(:,2);
         dN_deta.'*Node_lc(:,1) dN_deta.'*Node_lc(:,2)];
    
    dNi = J\[dN_dxi dN_deta].';
    
    Bm = zeros(3,8);
    
    Bm(1,1:2:7) = dNi(1,:);
    Bm(2,2:2:8) = dNi(2,:);
    Bm(3,1:2:7) = dNi(2,:);
    Bm(3,2:2:8) = dNi(1,:);
    
     Nex = [1-xi^2 -2*xi 0; %N5
            1-eta^2 0 -2*eta]; %N6
     dNex = J0\((det(J0)/det(J))*Nex(:,2:3).');
     Bex = zeros(3,4);
     Bex(1,[1 3])=dNex(1,:);
     Bex(2,[2 4])=dNex(2,:);
     Bex(3,[1 3])=dNex(2,:);
     Bex(3,[2 4])=dNex(1,:);
        
% KeLCS([1 2 7 8 13 14 19 20],[1 2 7 8 13 14 19 20]) = KeLCS([1 2 7 8 13 14 19 20],[1 2 7 8 13 14 19 20]) ...
%                                                      + Bm.'*Dm*Bm*Hm(1,iGauss)*det(J);
 
K1 = K1 + Bm.'*Dm*Bm*Hm(1,iGauss)*det(J);
K21 = K21 + Bex.'*Dm*Bm*Hm(1,iGauss)*det(J);
K2 = K2 + Bex.'*Dm*Bex*Hm(1,iGauss)*det(J);                                                
K12 = K12 + Bm.'*Dm*Bex*Hm(1,iGauss)*det(J);
                                                 
if nargout>1
    N = zeros(2,8);
    N(1,1:2:7) = Ni;
    N(2,2:2:8) = Ni;
MeLCS([1 2 7 8 13 14 19 20],[1 2 7 8 13 14 19 20]) = MeLCS([1 2 7 8 13 14 19 20],[1 2 7 8 13 14 19 20]) ...
                                                     + N.'*rho*h*N*Hm(1,iGauss)*det(J);
end

end

KeLCS([1 2 7 8 13 14 19 20],[1 2 7 8 13 14 19 20]) = K1-K12*(K2\K21);


% bending stiffness
% Node sequence of the four triangles
seq=[1 2 3;
     1 3 4;
     1 2 4;
     2 3 4];
if nargout > 1
for n =1:4
[Kt,Mt] = ke_dkt(Node_lc(seq(n,:),:),h,E,nu,rho);
    for indi =1:3
       for indj=1:3
KeLCS((3:5)+6*(seq(n,indi)-1),(3:5)+6*(seq(n,indj)-1)) = KeLCS((3:5)+6*(seq(n,indi)-1),(3:5)+6*(seq(n,indj)-1))...
                                                   + Kt((1:3)+3*(indi-1),(1:3)+3*(indj-1))/2;
MeLCS((3:5)+6*(seq(n,indi)-1),(3:5)+6*(seq(n,indj)-1)) = MeLCS((3:5)+6*(seq(n,indi)-1),(3:5)+6*(seq(n,indj)-1))...
                                                   + Mt((1:3)+3*(indi-1),(1:3)+3*(indj-1))/2;                                             
                                               
       end
    end
end
else
for n =1:4
Kt = ke_dkt(Node_lc(seq(n,:),:),h,E,nu);    
    for indi =1:3
       for indj=1:3
KeLCS((3:5)+6*(seq(n,indi)-1),(3:5)+6*(seq(n,indj)-1)) = KeLCS((3:5)+6*(seq(n,indi)-1),(3:5)+6*(seq(n,indj)-1))...
                                                   + Kt((1:3)+3*(indi-1),(1:3)+3*(indj-1))/2;
       end
    end
end
end

% eenvoudige oplossing voor drilling uit COOK p575
% Aelem = (1/2)*norm(cross((Node_lc(3,:)-Node_lc(1,:)),(Node_lc(4,:)-Node_lc(2,:))));
% KeLCS(6:6:24,6:6:24)=(1e-5)*E*Aelem*h(1)*[1 -1/3 -1/3 -1/3;
%                                           -1/3 1 -1/3 -1/3;
%                                           -1/3 -1/3 1 -1/3;
%                                           -1/3 -1/3 -1/3 1];

KeD = diag(KeLCS);
minKeD = min(KeD(KeD ~= 0));
for e=6:6:24
KeLCS(e,e) = minKeD/1000; % eenvoudige oplossing voor drilling uit BATHE p209                                      
end

%oplossing voor drilling volgens Kanok-Nukulchai
% kappa = 10;
% 
% Ni = 0.25*[1 1 1 1];
% J = [0.25*[-1 1 1 -1]*Node_lc(:,[1,2]);0.25*[-1 -1 1 1]*Node_lc(:,[1,2])];
% dNi = J\(0.25*[-1 1 1 -1; -1 -1 1 1]);
% 
% 
% Bt = zeros(1,24);
% Bt([1 7 13 19]) = 0.5*dNi(2,:);
% Bt([2 8 14 20]) = -0.5*dNi(1,:);
% Bt([6 12 18 24]) = Ni;
% 
% Kt = 4*kappa*(E/(2*(1+nu)))*h*(Bt.'*Bt)*det(J);
% 
% KeLCS = KeLCS+Kt;

end