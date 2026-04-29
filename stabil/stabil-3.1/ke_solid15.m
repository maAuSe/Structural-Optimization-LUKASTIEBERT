function [Ke,Me,dKedx]=ke_solid15(Node,Section,Material,Options,dNodedx,dSectiondx)

%KE_SOLID15   Compute the element stiffness and mass matrix for a solid15 element.
%
%   [Ke,Me] = ke_solid15(Node,Section,Material,Options)   computes element 
%   stiffness and mass matrix in the global coordinate system for a solid15 
%   prismatic element.
%
%   Node       Node definitions           [x y z] (15 * 3)
%              Nodes should have the following order:
%                    6
%                  / | \
%                12  |  11
%                /  15   \ 
%               /    |    \  
%              4 ---10 ----5
%              |     3     |
%              |   /   \   |   
%             13  9     8  14  
%              | /       \ | 
%              |/         \|
%              1 --- 7 --- 2
%
%   Section    Section definition         []
%   Material   Material definition        [E nu rho]
%   Options    Element options struct. Fields: []
%   Ke         Element stiffness matrix (45 * 45)
%   Me         Element mass matrix (45 * 45)

% Stijn Francois
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

% MATERIAL PROPERTIES
rho=Material(3);

C=cmat_isotropic('3d',Section,Material);

% --- Integration points and weights---
xiTri=[0.16666666666667  0.16666666666667
       0.16666666666667  0.66666666666667
       0.66666666666667  0.16666666666667];
HTri=[0.16666666666667
      0.16666666666667
      0.16666666666666];
nXiTri=size(xiTri,1);
xi1D=[0.774596669241484;0;-0.774596669241484];
H1D =[0.555555555555555;0.888888888888889;0.555555555555555];
nXi1D=numel(xi1D);


nXi=nXiTri*nXi1D;
xi=zeros(nXi,3);
H=zeros(nXi,1);
for iXi=1:nXiTri
  for jXi=1:nXi1D
    ind= nXi1D*(iXi-1)+jXi;
    xi(ind,1)= xiTri(iXi,1);            % xi1
    xi(ind,2)= xiTri(iXi,2);            % xi2
    xi(ind,3)= xi1D(jXi);               % xi3
    H(ind)=HTri(iXi)*H1D(jXi);
  end
end

% --- Shape function ---
N=zeros(nXi,15);
N(:, 1)=-(xi(:,3).^2 - 1).*(xi(:,1)/2 + xi(:,2)/2 - 1/2) - (xi(:,3) - 1).*(2.*xi(:,1) + 2.*xi(:,2) - 1).*(xi(:,1)/2 + xi(:,2)/2 - 1/2);
N(:, 2)=(xi(:,1).*(xi(:,3).^2 - 1))/2 - (xi(:,1).*(2.*xi(:,1) - 1).*(xi(:,3) - 1))/2;
N(:, 3)=(xi(:,2).*(xi(:,3).^2 - 1))/2 - (xi(:,2).*(2.*xi(:,2) - 1).*(xi(:,3) - 1))/2;
N(:, 4)=(xi(:,3) + 1).*(2.*xi(:,1) + 2.*xi(:,2) - 1).*(xi(:,1)/2 + xi(:,2)/2 - 1/2) - (xi(:,3).^2 - 1).*(xi(:,1)/2 + xi(:,2)/2 - 1/2);
N(:, 5)=(xi(:,1).*(xi(:,3).^2 - 1))/2 + (xi(:,1).*(2.*xi(:,1) - 1).*(xi(:,3) + 1))/2;
N(:, 6)=(xi(:,2).*(xi(:,3).^2 - 1))/2 + (xi(:,2).*(2.*xi(:,2) - 1).*(xi(:,3) + 1))/2;
N(:, 7)=xi(:,1).*(xi(:,3) - 1).*(2.*xi(:,1) + 2.*xi(:,2) - 2);
N(:, 8)=-2.*xi(:,1).*xi(:,2).*(xi(:,3) - 1);
N(:, 9)=2.*xi(:,2).*(xi(:,3) - 1).*(xi(:,1) + xi(:,2) - 1);
N(:,10)=-xi(:,1).*(xi(:,3) + 1).*(2.*xi(:,1) + 2.*xi(:,2) - 2);
N(:,11)=2.*xi(:,1).*xi(:,2).*(xi(:,3) + 1);
N(:,12)=-2.*xi(:,2).*(xi(:,3) + 1).*(xi(:,1) + xi(:,2) - 1);
N(:,13)=(xi(:,3).^2 - 1).*(xi(:,1) + xi(:,2) - 1);
N(:,14)=-xi(:,1).*(xi(:,3).^2 - 1);
N(:,15)=-xi(:,2).*(xi(:,3).^2 - 1);

% --- Shape function derivative ---
dN=zeros(nXi,45);
dN(:, 1)=-((xi(:,3) - 1).*(4.*xi(:,1) + 4.*xi(:,2) + xi(:,3) - 2))/2;                 % dN01(xi1,xi2,xi3)/dxi1
dN(:, 2)=((xi(:,3) - 1).*(xi(:,3) - 4.*xi(:,1) + 2))/2;                               % dN02(xi1,xi2,xi3)/dxi1
dN(:, 3)=0;                                                                           % dN03(xi1,xi2,xi3)/dxi1
dN(:, 4)=((xi(:,3) + 1).*(4.*xi(:,1) + 4.*xi(:,2) - xi(:,3) - 2))/2;                  % dN04(xi1,xi2,xi3)/dxi1
dN(:, 5)=((xi(:,3) + 1).*(4.*xi(:,1) + xi(:,3) - 2))/2;                               % dN05(xi1,xi2,xi3)/dxi1
dN(:, 6)=0;                                                                           % dN06(xi1,xi2,xi3)/dxi1
dN(:, 7)=2.*(xi(:,3) - 1).*(2.*xi(:,1) + xi(:,2) - 1);                                % dN07(xi1,xi2,xi3)/dxi1
dN(:, 8)=-2.*xi(:,2).*(xi(:,3) - 1);                                                  % dN08(xi1,xi2,xi3)/dxi1
dN(:, 9)=2.*xi(:,2).*(xi(:,3) - 1);                                                   % dN09(xi1,xi2,xi3)/dxi1
dN(:,10)=-2.*(xi(:,3) + 1).*(2.*xi(:,1) + xi(:,2) - 1);                               % dN10(xi1,xi2,xi3)/dxi1
dN(:,11)=2.*xi(:,2).*(xi(:,3) + 1);                                                   % dN11(xi1,xi2,xi3)/dxi1
dN(:,12)=-2.*xi(:,2).*(xi(:,3) + 1);                                                  % dN12(xi1,xi2,xi3)/dxi1
dN(:,13)=xi(:,3).^2 - 1;                                                              % dN13(xi1,xi2,xi3)/dxi1
dN(:,14)=1 - xi(:,3).^2;                                                              % dN14(xi1,xi2,xi3)/dxi1
dN(:,15)=0;                                                                           % dN15(xi1,xi2,xi3)/dxi1
dN(:,16)=-((xi(:,3) - 1).*(4.*xi(:,1) + 4.*xi(:,2) + xi(:,3) - 2))/2;                 % dN01(xi1,xi2,xi3)/dxi2
dN(:,17)=0;                                                                           % dN02(xi1,xi2,xi3)/dxi2
dN(:,18)=((xi(:,3) - 1).*(xi(:,3) - 4.*xi(:,2) + 2))/2;                               % dN03(xi1,xi2,xi3)/dxi2
dN(:,19)=((xi(:,3) + 1).*(4.*xi(:,1) + 4.*xi(:,2) - xi(:,3) - 2))/2;                  % dN04(xi1,xi2,xi3)/dxi2
dN(:,20)=0;                                                                           % dN05(xi1,xi2,xi3)/dxi2
dN(:,21)=((xi(:,3) + 1).*(4.*xi(:,2) + xi(:,3) - 2))/2;                               % dN06(xi1,xi2,xi3)/dxi2
dN(:,22)=2.*xi(:,1).*(xi(:,3) - 1);                                                   % dN07(xi1,xi2,xi3)/dxi2
dN(:,23)=-2.*xi(:,1).*(xi(:,3) - 1);                                                  % dN08(xi1,xi2,xi3)/dxi2
dN(:,24)=2.*(xi(:,3) - 1).*(xi(:,1) + 2.*xi(:,2) - 1);                                % dN09(xi1,xi2,xi3)/dxi2
dN(:,25)=-2.*xi(:,1).*(xi(:,3) + 1);                                                  % dN10(xi1,xi2,xi3)/dxi2
dN(:,26)=2.*xi(:,1).*(xi(:,3) + 1);                                                   % dN11(xi1,xi2,xi3)/dxi2
dN(:,27)=-2.*(xi(:,3) + 1).*(xi(:,1) + 2.*xi(:,2) - 1);                               % dN12(xi1,xi2,xi3)/dxi2
dN(:,28)=xi(:,3).^2 - 1;                                                              % dN13(xi1,xi2,xi3)/dxi2
dN(:,29)=0;                                                                           % dN14(xi1,xi2,xi3)/dxi2
dN(:,30)=1 - xi(:,3).^2;                                                              % dN15(xi1,xi2,xi3)/dxi2
dN(:,31)=-((xi(:,1) + xi(:,2) - 1).*(2.*xi(:,1) + 2.*xi(:,2) + 2.*xi(:,3) - 1))/2;    % dN01(xi1,xi2,xi3)/dxi3
dN(:,32)=(xi(:,1).*(2.*xi(:,3) - 2.*xi(:,1) + 1))/2;                                  % dN02(xi1,xi2,xi3)/dxi3
dN(:,33)=(xi(:,2).*(2.*xi(:,3) - 2.*xi(:,2) + 1))/2;                                  % dN03(xi1,xi2,xi3)/dxi3
dN(:,34)=((xi(:,1) + xi(:,2) - 1).*(2.*xi(:,1) + 2.*xi(:,2) - 2.*xi(:,3) - 1))/2;     % dN04(xi1,xi2,xi3)/dxi3
dN(:,35)=(xi(:,1).*(2.*xi(:,1) + 2.*xi(:,3) - 1))/2;                                  % dN05(xi1,xi2,xi3)/dxi3
dN(:,36)=(xi(:,2).*(2.*xi(:,2) + 2.*xi(:,3) - 1))/2;                                  % dN06(xi1,xi2,xi3)/dxi3
dN(:,37)=2.*xi(:,1).*(xi(:,1) + xi(:,2) - 1);                                         % dN07(xi1,xi2,xi3)/dxi3
dN(:,38)=-2.*xi(:,1).*xi(:,2);                                                        % dN08(xi1,xi2,xi3)/dxi3
dN(:,39)=2.*xi(:,2).*(xi(:,1) + xi(:,2) - 1);                                         % dN09(xi1,xi2,xi3)/dxi3
dN(:,40)=-2.*xi(:,1).*(xi(:,1) + xi(:,2) - 1);                                        % dN10(xi1,xi2,xi3)/dxi3
dN(:,41)=2.*xi(:,1).*xi(:,2);                                                         % dN11(xi1,xi2,xi3)/dxi3
dN(:,42)=-2.*xi(:,2).*(xi(:,1) + xi(:,2) - 1);                                        % dN12(xi1,xi2,xi3)/dxi3
dN(:,43)=2.*xi(:,3).*(xi(:,1) + xi(:,2) - 1);                                         % dN13(xi1,xi2,xi3)/dxi3
dN(:,44)=-2.*xi(:,1).*xi(:,3);                                                        % dN14(xi1,xi2,xi3)/dxi3
dN(:,45)=-2.*xi(:,2).*xi(:,3);                                                        % dN15(xi1,xi2,xi3)/dxi3

% --- Jacobian ---
% Jac =  [dx/dxi1 dx/dxi2 dx/dxi3
%         dy/dxi1 dy/dxi2 dy/dxi3
%         dz/dxi1 dz/dxi2 dz/dxi3]
Jac=zeros(nXi,3,3);
for iXi=1:nXi
    for iNod=1:15
      Jac(iXi,1,1)=Jac(iXi,1,1)+dN(iXi,iNod)*Node(iNod,1);
      Jac(iXi,1,2)=Jac(iXi,1,2)+dN(iXi,iNod)*Node(iNod,2);
      Jac(iXi,1,3)=Jac(iXi,1,3)+dN(iXi,iNod)*Node(iNod,3);
      Jac(iXi,2,1)=Jac(iXi,2,1)+dN(iXi,15+iNod)*Node(iNod,1);
      Jac(iXi,2,2)=Jac(iXi,2,2)+dN(iXi,15+iNod)*Node(iNod,2);
      Jac(iXi,2,3)=Jac(iXi,2,3)+dN(iXi,15+iNod)*Node(iNod,3);
      Jac(iXi,3,1)=Jac(iXi,3,1)+dN(iXi,30+iNod)*Node(iNod,1);
      Jac(iXi,3,2)=Jac(iXi,3,2)+dN(iXi,30+iNod)*Node(iNod,2);
      Jac(iXi,3,3)=Jac(iXi,3,3)+dN(iXi,30+iNod)*Node(iNod,3);
  end
end



% --- Integration: loop over integration points ---
Me=zeros(45,45);
Ke=zeros(45,45);
for iXi=1:nXi
  JacUtil=[Jac(iXi,1,1),Jac(iXi,1,2),Jac(iXi,1,3);Jac(iXi,2,1),Jac(iXi,2,2),Jac(iXi,2,3);Jac(iXi,3,1),Jac(iXi,3,2),Jac(iXi,3,3)];
  detJac=det(JacUtil);
  DN1= JacUtil\[dN(iXi, 1);dN(iXi,16);dN(iXi,31)];
  DN2= JacUtil\[dN(iXi, 2);dN(iXi,17);dN(iXi,32)];
  DN3= JacUtil\[dN(iXi, 3);dN(iXi,18);dN(iXi,33)];
  DN4= JacUtil\[dN(iXi, 4);dN(iXi,19);dN(iXi,34)];
  DN5= JacUtil\[dN(iXi, 5);dN(iXi,20);dN(iXi,35)];
  DN6= JacUtil\[dN(iXi, 6);dN(iXi,21);dN(iXi,36)];
  DN7= JacUtil\[dN(iXi, 7);dN(iXi,22);dN(iXi,37)];
  DN8= JacUtil\[dN(iXi, 8);dN(iXi,23);dN(iXi,38)];
  DN9= JacUtil\[dN(iXi, 9);dN(iXi,24);dN(iXi,39)];
  DN10=JacUtil\[dN(iXi,10);dN(iXi,25);dN(iXi,40)];
  DN11=JacUtil\[dN(iXi,11);dN(iXi,26);dN(iXi,41)];
  DN12=JacUtil\[dN(iXi,12);dN(iXi,27);dN(iXi,42)];
  DN13=JacUtil\[dN(iXi,13);dN(iXi,28);dN(iXi,43)];
  DN14=JacUtil\[dN(iXi,14);dN(iXi,29);dN(iXi,44)];
  DN15=JacUtil\[dN(iXi,15);dN(iXi,30);dN(iXi,45)];

  B=[DN1(1)  0        0      DN2(1)  0        0      DN3(1)  0        0      DN4(1)  0        0      DN5(1)  0        0      DN6(1)  0        0      DN7(1)  0        0      DN8(1)  0        0       DN9(1)  0        0      DN10(1)  0         0       DN11(1)  0         0       DN12(1)  0         0       DN13(1)  0         0       DN14(1)  0         0       DN15(1)  0         0
     0       DN1(2)   0      0       DN2(2)   0      0       DN3(2)   0      0       DN4(2)   0      0       DN5(2)   0      0       DN6(2)   0      0       DN7(2)   0      0       DN8(2)   0       0       DN9(2)   0      0        DN10(2)   0       0        DN11(2)   0       0        DN12(2)   0       0        DN13(2)   0       0        DN14(2)   0       0        DN15(2)   0
     0       0        DN1(3) 0       0        DN2(3) 0       0        DN3(3) 0       0        DN4(3) 0       0        DN5(3) 0       0        DN6(3) 0       0        DN7(3) 0       0        DN8(3)  0       0        DN9(3) 0        0         DN10(3) 0        0         DN11(3) 0        0         DN12(3) 0        0         DN13(3) 0        0         DN14(3) 0        0         DN15(3)
     DN1(2)  DN1(1)   0      DN2(2)  DN2(1)   0      DN3(2)  DN3(1)   0      DN4(2)  DN4(1)   0      DN5(2)  DN5(1)   0      DN6(2)  DN6(1)   0      DN7(2)  DN7(1)   0      DN8(2)  DN8(1)   0       DN9(2)  DN9(1)   0      DN10(2)  DN10(1)   0       DN11(2)  DN11(1)   0       DN12(2)  DN12(1)   0       DN13(2)  DN13(1)   0       DN14(2)  DN14(1)   0       DN15(2)  DN15(1)   0
     0       DN1(3)   DN1(2) 0       DN2(3)   DN2(2) 0       DN3(3)   DN3(2) 0       DN4(3)   DN4(2) 0       DN5(3)   DN5(2) 0       DN6(3)   DN6(2) 0       DN7(3)   DN7(2) 0       DN8(3)   DN8(2)  0       DN9(3)   DN9(2) 0        DN10(3)   DN10(2) 0        DN11(3)   DN11(2) 0        DN12(3)   DN12(2) 0        DN13(3)   DN13(2) 0        DN14(3)   DN14(2) 0        DN15(3)   DN15(2)
     DN1(3)  0        DN1(1) DN2(3)  0        DN2(1) DN3(3)  0        DN3(1) DN4(3)  0        DN4(1) DN5(3)  0        DN5(1) DN6(3)  0        DN6(1) DN7(3)  0        DN7(1) DN8(3)  0        DN8(1)  DN9(3)  0        DN9(1) DN10(3)  0         DN10(1) DN11(3)  0         DN11(1) DN12(3)  0         DN12(1) DN13(3)  0         DN13(1) DN14(3)  0         DN14(1) DN15(3)  0         DN15(1)];

  Nxi=[N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        0        N(iXi,9) 0        0        N(iXi,10) 0         0         N(iXi,11) 0        0        N(iXi,12) 0         0          N(iXi,13) 0        0          N(iXi,14) 0        0          N(iXi,15) 0        0
       0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        0        N(iXi,9) 0        0         N(iXi,10) 0         0        N(iXi,11) 0        0         N(iXi,12) 0          0        N(iXi,13) 0          0        N(iXi,14) 0          0        N(iXi,15) 0
       0        0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        0        N(iXi,9) 0         0         N(iXi,10) 0        0        N(iXi,11) 0         0         N(iXi,12)  0        0         N(iXi,13)  0        0         N(iXi,14)  0        0         N(iXi,15)];
  Ke=Ke+H(iXi)*detJac*(B.'*C*B);
  Me=Me+H(iXi)*detJac*(rho*Nxi.'*Nxi);
end
