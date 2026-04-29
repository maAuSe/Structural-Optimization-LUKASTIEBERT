function [Ke,Me,dKedx]=ke_solid20(Node,Section,Material,Options,dNodedx,dSectiondx)
%KE_SOLID20   Compute the element stiffness and mass matrix for a solid20 element.
%
%   [Ke,Me] = ke_solid20(Node,Section,Material,Options)   computes element 
%   stiffness and mass matrix in the global coordinate system for a solid20 
%   element.
%
%   Node       Node definitions           [x y z] (20 * 3)
%              Nodes should have the following order:
%                   8---15----7
%                  /|        /|
%                16 |      14 |
%                / 20      / 19 
%               /   |     /   |
%              5---13----6    |
%              |    4--11|----3
%              |   /     |   /
%             17 12     18  10
%              | /       | /
%              |/        |/
%              1----9----2
%
%   Section    Section definition         []
%   Material   Material definition        [E nu rho]
%   Options    Element options struct. Fields: []
%   Ke         Element stiffness matrix (60 * 60)
%   Me         Element mass matrix (60 * 60)
%
%   See also KE_BEAM, ASMKM, KE_TRUSS.

% Stijn Francois, Miche Jansen
% 2013

if nargin<5, dNodedx = []; end
if nargin<6, dSectiondx = []; end
if (nnz(dNodedx)+nnz(dSectiondx)) > 0
    error('Sensitivities have not been implemented yet.')
end
nVar = 0;
if nargout>2 
    if ~isempty(dNodedx) || ~isempty(dSectiondx)
        nVar = max(size(dNodedx,3),size(dSectiondx,3));
    end
    dKedx = cell(nVar,1);
end

nNode=size(unique(Node,'rows'),1);
if nNode==15  % Degenerate Prismatic Element
  [Ke,Me]=ke_solid15(Node([1 2 3 5 6 7 9 10 12 13 14 16 17 18 19],:),Section,Material,Options);
else

  % MATERIAL PROPERTIES
  if numel(Material)>2
    rho=Material(3);
  else
    rho=0;
  end
  
  % INITIALIZATION
  Me=zeros(60,60);
  Ke=zeros(60,60);
  
  C=cmat_isotropic('3d',Section,Material);
  
  % --- Integration points and weights---
  xi1D=[0.774596669241484;0;-0.774596669241484];
  H1D =[0.555555555555555;0.888888888888889;0.555555555555555];
  nXi1D=numel(xi1D);
  
  nXi=nXi1D^3;
  xi=zeros(nXi,3);
  H=zeros(nXi,1);
  for iXi=1:nXi1D
    for jXi=1:nXi1D
      for kXi=1:nXi1D
        ind=nXi1D^2*(iXi-1)+nXi1D*(jXi-1)+kXi;
        xi(ind,1)= xi1D(iXi);              % xi1
        xi(ind,2)= xi1D(jXi);              % xi2
        xi(ind,3)= xi1D(kXi);              % xi3
        H(ind)=H1D(iXi)*H1D(jXi)*H1D(kXi);
      end
    end
  end
  
  % --- Shape function ---
  N=zeros(nXi,20);
  N(:, 1)= 0.125.*((xi(:,1)-1).*(xi(:,2)-1).*(xi(:,3)-1).*(xi(:,1)+xi(:,2)+xi(:,3)+2));  % N01(xi1,xi2,xi3)
  N(:, 2)=-0.125.*((xi(:,1)+1).*(xi(:,2)-1).*(xi(:,3)-1).*(xi(:,2)-xi(:,1)+xi(:,3)+2));  % N02(xi1,xi2,xi3)
  N(:, 3)=-0.125.*((xi(:,1)+1).*(xi(:,2)+1).*(xi(:,3)-1).*(xi(:,1)+xi(:,2)-xi(:,3)-2));  % N03(xi1,xi2,xi3)
  N(:, 4)=-0.125.*((xi(:,1)-1).*(xi(:,2)+1).*(xi(:,3)-1).*(xi(:,1)-xi(:,2)+xi(:,3)+2));  % N04(xi1,xi2,xi3)
  N(:, 5)=-0.125.*((xi(:,1)-1).*(xi(:,2)-1).*(xi(:,3)+1).*(xi(:,1)+xi(:,2)-xi(:,3)+2));  % N05(xi1,xi2,xi3)
  N(:, 6)=-0.125.*((xi(:,1)+1).*(xi(:,2)-1).*(xi(:,3)+1).*(xi(:,1)-xi(:,2)+xi(:,3)-2));  % N06(xi1,xi2,xi3)
  N(:, 7)= 0.125.*((xi(:,1)+1).*(xi(:,2)+1).*(xi(:,3)+1).*(xi(:,1)+xi(:,2)+xi(:,3)-2));  % N07(xi1,xi2,xi3)
  N(:, 8)= 0.125.*((xi(:,1)-1).*(xi(:,2)+1).*(xi(:,3)+1).*(xi(:,1)-xi(:,2)-xi(:,3)+2));  % N08(xi1,xi2,xi3)
  N(:, 9)=-0.25.*(xi(:,1).^2 - 1).*(xi(:,2) - 1)   .*(xi(:,3) - 1);                      % N09(xi1,xi2,xi3)
  N(:,10)= 0.25.*(xi(:,1) + 1)   .*(xi(:,2).^2 - 1).*(xi(:,3) - 1);                      % N10(xi1,xi2,xi3)
  N(:,11)= 0.25.*(xi(:,1).^2 - 1).*(xi(:,2) + 1)   .*(xi(:,3) - 1);                      % N11(xi1,xi2,xi3)
  N(:,12)=-0.25.*(xi(:,1) - 1)   .*(xi(:,2).^2 - 1).*(xi(:,3) - 1);                      % N12(xi1,xi2,xi3)
  N(:,13)= 0.25.*(xi(:,1).^2 - 1).*(xi(:,2) - 1)   .*(xi(:,3) + 1);                      % N13(xi1,xi2,xi3)
  N(:,14)=-0.25.*(xi(:,1) + 1)   .*(xi(:,2).^2 - 1).*(xi(:,3) + 1);                      % N14(xi1,xi2,xi3)
  N(:,15)=-0.25.*(xi(:,1).^2 - 1).*(xi(:,2) + 1)   .*(xi(:,3) + 1);                      % N15(xi1,xi2,xi3)
  N(:,16)= 0.25.*(xi(:,1) - 1)   .*(xi(:,2).^2 - 1).*(xi(:,3) + 1);                      % N16(xi1,xi2,xi3)
  N(:,17)=-0.25.*(xi(:,1) - 1)   .*(xi(:,3).^2 - 1).*(xi(:,2) - 1);                      % N17(xi1,xi2,xi3)
  N(:,18)= 0.25.*(xi(:,1) + 1)   .*(xi(:,3).^2 - 1).*(xi(:,2) - 1);                      % N18(xi1,xi2,xi3)
  N(:,19)=-0.25.*(xi(:,1) + 1)   .*(xi(:,3).^2 - 1).*(xi(:,2) + 1);                      % N19(xi1,xi2,xi3)
  N(:,20)= 0.25.*(xi(:,1) - 1)   .*(xi(:,3).^2 - 1).*(xi(:,2) + 1);                      % N20(xi1,xi2,xi3)
  
  
  % --- Shape function derivative ---
  dN=zeros(nXi,60);
  dN(:, 1)= ((xi(:,2) - 1).*(xi(:,3) - 1).*(2.*xi(:,1) + xi(:,2) + xi(:,3) + 1))/8;  % dN01(xi1,xi2,xi3)/dxi1
  dN(:, 2)=-((xi(:,2) - 1).*(xi(:,3) - 1).*(xi(:,2) - 2.*xi(:,1) + xi(:,3) + 1))/8;  % dN02(xi1,xi2,xi3)/dxi1
  dN(:, 3)=-((xi(:,2) + 1).*(xi(:,3) - 1).*(2.*xi(:,1) + xi(:,2) - xi(:,3) - 1))/8;  % dN03(xi1,xi2,xi3)/dxi1
  dN(:, 4)=-((xi(:,2) + 1).*(xi(:,3) - 1).*(2.*xi(:,1) - xi(:,2) + xi(:,3) + 1))/8;  % dN04(xi1,xi2,xi3)/dxi1
  dN(:, 5)=-((xi(:,2) - 1).*(xi(:,3) + 1).*(2.*xi(:,1) + xi(:,2) - xi(:,3) + 1))/8;  % dN05(xi1,xi2,xi3)/dxi1
  dN(:, 6)=-((xi(:,2) - 1).*(xi(:,3) + 1).*(2.*xi(:,1) - xi(:,2) + xi(:,3) - 1))/8;  % dN06(xi1,xi2,xi3)/dxi1
  dN(:, 7)= ((xi(:,2) + 1).*(xi(:,3) + 1).*(2.*xi(:,1) + xi(:,2) + xi(:,3) - 1))/8;  % dN07(xi1,xi2,xi3)/dxi1
  dN(:, 8)= ((xi(:,2) + 1).*(xi(:,3) + 1).*(2.*xi(:,1) - xi(:,2) - xi(:,3) + 1))/8;  % dN08(xi1,xi2,xi3)/dxi1
  dN(:, 9)= -(xi(:,1).*(xi(:,2) - 1).*(xi(:,3) - 1))/2;                              % dN09(xi1,xi2,xi3)/dxi1
  dN(:,10)= ((xi(:,2).^2 - 1).*(xi(:,3) - 1))/4;                                     % dN10(xi1,xi2,xi3)/dxi1
  dN(:,11)=  (xi(:,1).*(xi(:,2) + 1).*(xi(:,3) - 1))/2;                              % dN11(xi1,xi2,xi3)/dxi1
  dN(:,12)=-((xi(:,2).^2 - 1).*(xi(:,3) - 1))/4;                                     % dN12(xi1,xi2,xi3)/dxi1
  dN(:,13)=  (xi(:,1).*(xi(:,2) - 1).*(xi(:,3) + 1))/2;                              % dN13(xi1,xi2,xi3)/dxi1
  dN(:,14)=-((xi(:,2).^2 - 1).*(xi(:,3) + 1))/4;                                     % dN14(xi1,xi2,xi3)/dxi1
  dN(:,15)=-(xi(:,1).*(xi(:,2) + 1).*(xi(:,3) + 1))/2;                               % dN15(xi1,xi2,xi3)/dxi1
  dN(:,16)=((xi(:,2).^2 - 1).*(xi(:,3) + 1))/4;                                      % dN16(xi1,xi2,xi3)/dxi1
  dN(:,17)=-((xi(:,3).^2 - 1).*(xi(:,2) - 1))/4;                                     % dN17(xi1,xi2,xi3)/dxi1
  dN(:,18)=((xi(:,3).^2 - 1).*(xi(:,2) - 1))/4;                                      % dN18(xi1,xi2,xi3)/dxi1
  dN(:,19)=-((xi(:,3).^2 - 1).*(xi(:,2) + 1))/4;                                     % dN19(xi1,xi2,xi3)/dxi1
  dN(:,20)=((xi(:,3).^2 - 1).*(xi(:,2) + 1))/4;                                      % dN20(xi1,xi2,xi3)/dxi1
  dN(:,21)=((xi(:,1) - 1).*(xi(:,3) - 1).*(xi(:,1) + 2.*xi(:,2) + xi(:,3) + 1))/8;   % dN01(xi1,xi2,xi3)/dxi2
  dN(:,22)=-((xi(:,1) + 1).*(xi(:,3) - 1).*(2.*xi(:,2) - xi(:,1) + xi(:,3) + 1))/8;  % dN02(xi1,xi2,xi3)/dxi2
  dN(:,23)=-((xi(:,1) + 1).*(xi(:,3) - 1).*(xi(:,1) + 2.*xi(:,2) - xi(:,3) - 1))/8;  % dN03(xi1,xi2,xi3)/dxi2
  dN(:,24)=-((xi(:,1) - 1).*(xi(:,3) - 1).*(xi(:,1) - 2.*xi(:,2) + xi(:,3) + 1))/8;  % dN04(xi1,xi2,xi3)/dxi2
  dN(:,25)=-((xi(:,1) - 1).*(xi(:,3) + 1).*(xi(:,1) + 2.*xi(:,2) - xi(:,3) + 1))/8;  % dN05(xi1,xi2,xi3)/dxi2
  dN(:,26)=-((xi(:,1) + 1).*(xi(:,3) + 1).*(xi(:,1) - 2.*xi(:,2) + xi(:,3) - 1))/8;  % dN06(xi1,xi2,xi3)/dxi2
  dN(:,27)=((xi(:,1) + 1).*(xi(:,3) + 1).*(xi(:,1) + 2.*xi(:,2) + xi(:,3) - 1))/8;   % dN07(xi1,xi2,xi3)/dxi2
  dN(:,28)=((xi(:,1) - 1).*(xi(:,3) + 1).*(xi(:,1) - 2.*xi(:,2) - xi(:,3) + 1))/8;   % dN08(xi1,xi2,xi3)/dxi2
  dN(:,29)=-((xi(:,1).^2 - 1).*(xi(:,3) - 1))/4;                                     % dN09(xi1,xi2,xi3)/dxi2
  dN(:,30)=(xi(:,2).*(xi(:,1) + 1).*(xi(:,3) - 1))/2;                                % dN10(xi1,xi2,xi3)/dxi2
  dN(:,31)=((xi(:,1).^2 - 1).*(xi(:,3) - 1))/4;                                      % dN11(xi1,xi2,xi3)/dxi2
  dN(:,32)=-(xi(:,2).*(xi(:,1) - 1).*(xi(:,3) - 1))/2;                               % dN12(xi1,xi2,xi3)/dxi2
  dN(:,33)=((xi(:,1).^2 - 1).*(xi(:,3) + 1))/4;                                      % dN13(xi1,xi2,xi3)/dxi2
  dN(:,34)=-(xi(:,2).*(xi(:,1) + 1).*(xi(:,3) + 1))/2;                               % dN14(xi1,xi2,xi3)/dxi2
  dN(:,35)=-((xi(:,1).^2 - 1).*(xi(:,3) + 1))/4;                                     % dN15(xi1,xi2,xi3)/dxi2
  dN(:,36)=(xi(:,2).*(xi(:,1) - 1).*(xi(:,3) + 1))/2;                                % dN16(xi1,xi2,xi3)/dxi2
  dN(:,37)=-((xi(:,3).^2 - 1).*(xi(:,1) - 1))/4;                                     % dN17(xi1,xi2,xi3)/dxi2
  dN(:,38)=((xi(:,3).^2 - 1).*(xi(:,1) + 1))/4;                                      % dN18(xi1,xi2,xi3)/dxi2
  dN(:,39)=-((xi(:,3).^2 - 1).*(xi(:,1) + 1))/4;                                     % dN19(xi1,xi2,xi3)/dxi2
  dN(:,40)=((xi(:,3).^2 - 1).*(xi(:,1) - 1))/4;                                      % dN20(xi1,xi2,xi3)/dxi2
  dN(:,41)=((xi(:,1) - 1).*(xi(:,2) - 1).*(xi(:,1) + xi(:,2) + 2.*xi(:,3) + 1))/8;   % dN01(xi1,xi2,xi3)/dxi3
  dN(:,42)=-((xi(:,1) + 1).*(xi(:,2) - 1).*(xi(:,2) - xi(:,1) + 2.*xi(:,3) + 1))/8;  % dN02(xi1,xi2,xi3)/dxi3
  dN(:,43)=-((xi(:,1) + 1).*(xi(:,2) + 1).*(xi(:,1) + xi(:,2) - 2.*xi(:,3) - 1))/8;  % dN03(xi1,xi2,xi3)/dxi3
  dN(:,44)=-((xi(:,1) - 1).*(xi(:,2) + 1).*(xi(:,1) - xi(:,2) + 2.*xi(:,3) + 1))/8;  % dN04(xi1,xi2,xi3)/dxi3
  dN(:,45)=-((xi(:,1) - 1).*(xi(:,2) - 1).*(xi(:,1) + xi(:,2) - 2.*xi(:,3) + 1))/8;  % dN05(xi1,xi2,xi3)/dxi3
  dN(:,46)=-((xi(:,1) + 1).*(xi(:,2) - 1).*(xi(:,1) - xi(:,2) + 2.*xi(:,3) - 1))/8;  % dN06(xi1,xi2,xi3)/dxi3
  dN(:,47)=((xi(:,1) + 1).*(xi(:,2) + 1).*(xi(:,1) + xi(:,2) + 2.*xi(:,3) - 1))/8;   % dN07(xi1,xi2,xi3)/dxi3
  dN(:,48)=((xi(:,1) - 1).*(xi(:,2) + 1).*(xi(:,1) - xi(:,2) - 2.*xi(:,3) + 1))/8;   % dN08(xi1,xi2,xi3)/dxi3
  dN(:,49)=-((xi(:,1).^2 - 1).*(xi(:,2) - 1))/4;                                     % dN09(xi1,xi2,xi3)/dxi3
  dN(:,50)=((xi(:,2).^2 - 1).*(xi(:,1) + 1))/4;                                      % dN10(xi1,xi2,xi3)/dxi3
  dN(:,51)=((xi(:,1).^2 - 1).*(xi(:,2) + 1))/4;                                      % dN11(xi1,xi2,xi3)/dxi3
  dN(:,52)=-((xi(:,2).^2 - 1).*(xi(:,1) - 1))/4;                                     % dN12(xi1,xi2,xi3)/dxi3
  dN(:,53)=((xi(:,1).^2 - 1).*(xi(:,2) - 1))/4;                                      % dN13(xi1,xi2,xi3)/dxi3
  dN(:,54)=-((xi(:,2).^2 - 1).*(xi(:,1) + 1))/4;                                     % dN14(xi1,xi2,xi3)/dxi3
  dN(:,55)=-((xi(:,1).^2 - 1).*(xi(:,2) + 1))/4;                                     % dN15(xi1,xi2,xi3)/dxi3
  dN(:,56)=((xi(:,2).^2 - 1).*(xi(:,1) - 1))/4;                                      % dN16(xi1,xi2,xi3)/dxi3
  dN(:,57)=-(xi(:,3).*(xi(:,1) - 1).*(xi(:,2) - 1))/2;                               % dN17(xi1,xi2,xi3)/dxi3
  dN(:,58)=(xi(:,3).*(xi(:,1) + 1).*(xi(:,2) - 1))/2;                                % dN18(xi1,xi2,xi3)/dxi3
  dN(:,59)=-(xi(:,3).*(xi(:,1) + 1).*(xi(:,2) + 1))/2;                               % dN19(xi1,xi2,xi3)/dxi3
  dN(:,60)=(xi(:,3).*(xi(:,1) - 1).*(xi(:,2) + 1))/2;                                % dN20(xi1,xi2,xi3)/dxi3
  
  
  % --- Jacobian ---
  % Jac =  [dx/dxi1 dx/dxi2 dx/dxi3
  %         dy/dxi1 dy/dxi2 dy/dxi3
  %         dz/dxi1 dz/dxi2 dz/dxi3]
  Jac=zeros(nXi,3,3);
  for iXi=1:nXi
  	for iNod=1:20
  	  Jac(iXi,1,1)=Jac(iXi,1,1)+dN(iXi,iNod)*Node(iNod,1);  
  	  Jac(iXi,1,2)=Jac(iXi,1,2)+dN(iXi,iNod)*Node(iNod,2);  
      Jac(iXi,1,3)=Jac(iXi,1,3)+dN(iXi,iNod)*Node(iNod,3);  
  	  Jac(iXi,2,1)=Jac(iXi,2,1)+dN(iXi,20+iNod)*Node(iNod,1);
  	  Jac(iXi,2,2)=Jac(iXi,2,2)+dN(iXi,20+iNod)*Node(iNod,2);
      Jac(iXi,2,3)=Jac(iXi,2,3)+dN(iXi,20+iNod)*Node(iNod,3);
      Jac(iXi,3,1)=Jac(iXi,3,1)+dN(iXi,40+iNod)*Node(iNod,1);
  	  Jac(iXi,3,2)=Jac(iXi,3,2)+dN(iXi,40+iNod)*Node(iNod,2);
      Jac(iXi,3,3)=Jac(iXi,3,3)+dN(iXi,40+iNod)*Node(iNod,3);
    end
  end
  
  % --- Integration: loop over integration points ---
  for iXi=1:nXi
    JacUtil=[Jac(iXi,1,1),Jac(iXi,1,2),Jac(iXi,1,3);Jac(iXi,2,1),Jac(iXi,2,2),Jac(iXi,2,3);Jac(iXi,3,1),Jac(iXi,3,2),Jac(iXi,3,3)];
    detJac=det(JacUtil);
    DN1= JacUtil\[dN(iXi, 1);dN(iXi,21);dN(iXi,41)];
    DN2= JacUtil\[dN(iXi, 2);dN(iXi,22);dN(iXi,42)];
    DN3= JacUtil\[dN(iXi, 3);dN(iXi,23);dN(iXi,43)];
    DN4= JacUtil\[dN(iXi, 4);dN(iXi,24);dN(iXi,44)];
    DN5= JacUtil\[dN(iXi, 5);dN(iXi,25);dN(iXi,45)];
    DN6= JacUtil\[dN(iXi, 6);dN(iXi,26);dN(iXi,46)];
    DN7= JacUtil\[dN(iXi, 7);dN(iXi,27);dN(iXi,47)];
    DN8= JacUtil\[dN(iXi, 8);dN(iXi,28);dN(iXi,48)];
    DN9= JacUtil\[dN(iXi, 9);dN(iXi,29);dN(iXi,49)];
    DN10=JacUtil\[dN(iXi,10);dN(iXi,30);dN(iXi,50)];
    DN11=JacUtil\[dN(iXi,11);dN(iXi,31);dN(iXi,51)];
    DN12=JacUtil\[dN(iXi,12);dN(iXi,32);dN(iXi,52)];
    DN13=JacUtil\[dN(iXi,13);dN(iXi,33);dN(iXi,53)];
    DN14=JacUtil\[dN(iXi,14);dN(iXi,34);dN(iXi,54)];
    DN15=JacUtil\[dN(iXi,15);dN(iXi,35);dN(iXi,55)];
    DN16=JacUtil\[dN(iXi,16);dN(iXi,36);dN(iXi,56)];
    DN17=JacUtil\[dN(iXi,17);dN(iXi,37);dN(iXi,57)];
    DN18=JacUtil\[dN(iXi,18);dN(iXi,38);dN(iXi,58)];
    DN19=JacUtil\[dN(iXi,19);dN(iXi,39);dN(iXi,59)];
    DN20=JacUtil\[dN(iXi,20);dN(iXi,40);dN(iXi,60)];
  
    B=[DN1(1)  0        0      DN2(1)  0        0      DN3(1)  0        0      DN4(1)  0        0      DN5(1)  0        0      DN6(1)  0        0      DN7(1)  0        0      DN8(1)  0        0       DN9(1)  0        0      DN10(1)  0         0       DN11(1)  0         0       DN12(1)  0         0       DN13(1)  0         0       DN14(1)  0         0       DN15(1)  0         0       DN16(1)  0         0       DN17(1)  0         0       DN18(1)  0         0       DN19(1)  0         0       DN20(1)  0         0      
       0       DN1(2)   0      0       DN2(2)   0      0       DN3(2)   0      0       DN4(2)   0      0       DN5(2)   0      0       DN6(2)   0      0       DN7(2)   0      0       DN8(2)   0       0       DN9(2)   0      0        DN10(2)   0       0        DN11(2)   0       0        DN12(2)   0       0        DN13(2)   0       0        DN14(2)   0       0        DN15(2)   0       0        DN16(2)   0       0        DN17(2)   0       0        DN18(2)   0       0        DN19(2)   0       0        DN20(2)   0      
       0       0        DN1(3) 0       0        DN2(3) 0       0        DN3(3) 0       0        DN4(3) 0       0        DN5(3) 0       0        DN6(3) 0       0        DN7(3) 0       0        DN8(3)  0       0        DN9(3) 0        0         DN10(3) 0        0         DN11(3) 0        0         DN12(3) 0        0         DN13(3) 0        0         DN14(3) 0        0         DN15(3) 0        0         DN16(3) 0        0         DN17(3) 0        0         DN18(3) 0        0         DN19(3) 0        0         DN20(3)
       DN1(2)  DN1(1)   0      DN2(2)  DN2(1)   0      DN3(2)  DN3(1)   0      DN4(2)  DN4(1)   0      DN5(2)  DN5(1)   0      DN6(2)  DN6(1)   0      DN7(2)  DN7(1)   0      DN8(2)  DN8(1)   0       DN9(2)  DN9(1)   0      DN10(2)  DN10(1)   0       DN11(2)  DN11(1)   0       DN12(2)  DN12(1)   0       DN13(2)  DN13(1)   0       DN14(2)  DN14(1)   0       DN15(2)  DN15(1)   0       DN16(2)  DN16(1)   0       DN17(2)  DN17(1)   0       DN18(2)  DN18(1)   0       DN19(2)  DN19(1)   0       DN20(2)  DN20(1)   0      
       0       DN1(3)   DN1(2) 0       DN2(3)   DN2(2) 0       DN3(3)   DN3(2) 0       DN4(3)   DN4(2) 0       DN5(3)   DN5(2) 0       DN6(3)   DN6(2) 0       DN7(3)   DN7(2) 0       DN8(3)   DN8(2)  0       DN9(3)   DN9(2) 0        DN10(3)   DN10(2) 0        DN11(3)   DN11(2) 0        DN12(3)   DN12(2) 0        DN13(3)   DN13(2) 0        DN14(3)   DN14(2) 0        DN15(3)   DN15(2) 0        DN16(3)   DN16(2) 0        DN17(3)   DN17(2) 0        DN18(3)   DN18(2) 0        DN19(3)   DN19(2) 0        DN20(3)   DN20(2)
       DN1(3)  0        DN1(1) DN2(3)  0        DN2(1) DN3(3)  0        DN3(1) DN4(3)  0        DN4(1) DN5(3)  0        DN5(1) DN6(3)  0        DN6(1) DN7(3)  0        DN7(1) DN8(3)  0        DN8(1)  DN9(3)  0        DN9(1) DN10(3)  0         DN10(1) DN11(3)  0         DN11(1) DN12(3)  0         DN12(1) DN13(3)  0         DN13(1) DN14(3)  0         DN14(1) DN15(3)  0         DN15(1) DN16(3)  0         DN16(1) DN17(3)  0         DN17(1) DN18(3)  0         DN18(1) DN19(3)  0         DN19(1) DN20(3)  0         DN20(1)];
  
    Nxi=[N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        0        N(iXi,9) 0        0        N(iXi,10) 0         0         N(iXi,11) 0        0        N(iXi,12) 0         0          N(iXi,13) 0        0          N(iXi,14) 0        0          N(iXi,15) 0        0          N(iXi,16) 0        0          N(iXi,17) 0         0          N(iXi,18) 0         0          N(iXi,19) 0         0          N(iXi,20) 0         0        
         0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        0        N(iXi,9) 0        0         N(iXi,10) 0         0        N(iXi,11) 0        0         N(iXi,12) 0          0        N(iXi,13) 0          0        N(iXi,14) 0          0        N(iXi,15) 0          0        N(iXi,16) 0          0         N(iXi,17) 0          0         N(iXi,18) 0          0         N(iXi,19) 0          0         N(iXi,20) 0        
         0        0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        0        N(iXi,9) 0         0         N(iXi,10) 0        0        N(iXi,11) 0         0         N(iXi,12)  0        0         N(iXi,13)  0        0         N(iXi,14)  0        0         N(iXi,15)  0        0         N(iXi,16)  0         0         N(iXi,17)  0         0         N(iXi,18)  0         0         N(iXi,19)  0         0         N(iXi,20)];
    Ke=Ke+H(iXi)*detJac*(B.'*C*B);
    Me=Me+H(iXi)*detJac*(rho*Nxi.'*Nxi);
  end
end