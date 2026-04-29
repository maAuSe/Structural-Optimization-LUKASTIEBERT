function SeLCS = selcs_shell4(Node_lc,h,E,nu,UeLCS,Options)

%SELCS_SHELL4   Compute the element stresses for a shell4 element.
%
%   [SeLCS] = selcs_shell4(Node,Section,Material,UeGCS,Options)
%   computes the element stresses in the global and the  
%   local coordinate system for the shell4 element.
%
%   Node       Node definitions           [x y z] (4 * 3)
%   h          Shell thickness         
%   E          Young's modulus
%   nu         Poisson coefficient
%   rho        Mass density
%   UeLCS      Displacements (24 * nTimeSteps)
%   Options    Element options            {Option1 Option2 ...}
%   SeLCS      Element stresses in LCS in corner nodes IJKL and 
%              at top/mid/bot of shell (72 * nTimeSteps)   
%                                        [sxx syy szz sxy syz sxz]
%
%   See also ELEMSTRESS, SE_SHELL4.

% Miche Jansen
% 2010

%membrane
[xm,Hm]= gaussq(2);

D = E/(1-nu^2)*[1 nu 0;
     nu 1 0;
     0 0 (1-nu)/2];

SmeLCSg = zeros(12,size(UeLCS,2));
 

K21 = zeros(4,8);
K2 = zeros(4,4);


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
        
K21 = K21 + Bex.'*Dm*Bm*Hm(1,iGauss)*det(J);
K2 = K2 + Bex.'*Dm*Bex*Hm(1,iGauss)*det(J);                                                

end

Uex = -K2\K21*UeLCS([1 2 7 8 13 14 19 20],:);
 
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
    
    SmeLCSg((1:3)+3*(iGauss-1),:) = D*Bm*UeLCS([1 2 7 8 13 14 19 20],:)+D*Bex*Uex;
end
g=1/abs(xm(1,1))/2;
extrap = [(1+g)*eye(3) -0.5*eye(3) (1-g)*eye(3) -0.5*eye(3);...
          -0.5*eye(3) (1+g)*eye(3) -0.5*eye(3) (1-g)*eye(3);...
          (1-g)*eye(3) -0.5*eye(3) (1+g)*eye(3) -0.5*eye(3);...
          -0.5*eye(3) (1-g)*eye(3) -0.5*eye(3) (1+g)*eye(3)];
SmeLCS = extrap*SmeLCSg;

%bending
SbeLCS = zeros(40,size(UeLCS,2));

% Node sequence of the four triangles
seq=[1 2 3;
     3 4 1;
     4 1 2;
     2 3 4];
L = [1 0 0;1 1 0;1 0 1]; 

zeta = [1 -1];


for n =1:4 
Nodet = Node_lc(seq(n,:),:);
ind = [6*(seq(n,1)-1)+(3:5),6*(seq(n,2)-1)+(3:5),6*(seq(n,3)-1)+(3:5)];
UeLCSt = UeLCS(ind,:);
     
b = [Nodet(2,2)-Nodet(3,2);
     Nodet(3,2)-Nodet(1,2);
     Nodet(1,2)-Nodet(2,2)];
c = [Nodet(3,1)-Nodet(2,1);
     Nodet(1,1)-Nodet(3,1);
     Nodet(2,1)-Nodet(1,1)];

deter = b(2)*c(3)-b(3)*c(2);

Q = q_dkt(b,c,deter);

    
chi1 = L*Q((1:3),:)*UeLCSt/deter;
chi2 = L*Q((4:6),:)*UeLCSt/deter;
chi12 = L*Q((7:9),:)*UeLCSt/deter;

for z = 1:2 %itereren over de doorsnede (top,mid,bot)
    
%dwarskracht/h !
vx = h^2/12*(D(1,1)*b(2:3).'*Q(2:3,:)+D(1,2)*b(2:3).'*Q(5:6,:)+D(3,3)*c(2:3).'*Q(8:9,:))*UeLCSt/deter^2;
vy = h^2/12*(D(2,1)*c(2:3).'*Q(2:3,:)+D(2,2)*c(2:3).'*Q(5:6,:)+D(3,3)*b(2:3).'*Q(8:9,:))*UeLCSt/deter^2;


SbeLCS(5*(seq(n,2)-1)+(1:3)+20*(z-1),:) = (zeta(z)*h/2)*D*[chi1(2,:);chi2(2,:);chi12(2,:)];

SbeLCS(5*(seq(n,2)-1)+4+20*(z-1),:) = vy;
SbeLCS(5*(seq(n,2)-1)+5+20*(z-1),:) = vx; 
end
end

% SeLCS = zeros(72,size(UeLCS,2)); 
% SeLCS(1:6:19,:) = SmeLCS(1:3:10,:)+SbeLCS(1:5:16,:);
% SeLCS(25:6:43,:) = SmeLCS(1:3:10,:);
% SeLCS(49:6:67,:) = SmeLCS(1:3:10,:)+SbeLCS(21:5:36,:);
% SeLCS(2:6:20,:) = SmeLCS(2:3:11,:)+SbeLCS(2:5:17,:);
% SeLCS(26:6:44,:) = SmeLCS(2:3:11,:);
% SeLCS(50:6:68,:) = SmeLCS(2:3:11,:)+SbeLCS(22:5:37,:);
% SeLCS(4:6:22,:) = SmeLCS(3:3:12,:)+SbeLCS(3:5:18,:);
% SeLCS(28:6:46,:) = SmeLCS(3:3:12,:);
% SeLCS(52:6:70,:) = SmeLCS(3:3:12,:)+SbeLCS(23:5:38,:);
% SeLCS(5:6:23,:) = SbeLCS(4:5:19);
% SeLCS(29:6:47,:) = SbeLCS(24:5:39,:);
% SeLCS(53:6:71,:) = SbeLCS(24:5:39,:);
% SeLCS(6:6:24,:) = SbeLCS(5:5:20);
% SeLCS(30:6:48,:) = SbeLCS(25:5:40,:);
% SeLCS(54:6:72,:) = SbeLCS(25:5:40,:);

SeLCS = zeros(72,size(UeLCS,2)); 
SeLCS(1:6:19,:) = SmeLCS(1:3:10,:)+SbeLCS(1:5:16,:);
SeLCS(25:6:43,:) = SmeLCS(1:3:10,:);
SeLCS(49:6:67,:) = SmeLCS(1:3:10,:)+SbeLCS(21:5:36,:);
SeLCS(2:6:20,:) = SmeLCS(2:3:11,:)+SbeLCS(2:5:17,:);
SeLCS(26:6:44,:) = SmeLCS(2:3:11,:);
SeLCS(50:6:68,:) = SmeLCS(2:3:11,:)+SbeLCS(22:5:37,:);
SeLCS(4:6:22,:) = SmeLCS(3:3:12,:)+SbeLCS(3:5:18,:);
SeLCS(28:6:46,:) = SmeLCS(3:3:12,:);
SeLCS(52:6:70,:) = SmeLCS(3:3:12,:)+SbeLCS(23:5:38,:);
SeLCS(29:6:47,:) = 3/2*SbeLCS(4:5:19,:);
SeLCS(30:6:48,:) = 3/2*SbeLCS(5:5:20,:);

end


 
 




