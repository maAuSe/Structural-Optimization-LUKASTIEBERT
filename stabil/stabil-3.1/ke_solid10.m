function [Ke,Me,dKedx]=ke_solid10(Node,Section,Material,Options,dNodedx,dSectiondx)
%KE_SOLID10   Compute the element stiffness and mass matrix for a solid10 element.
%
%   [Ke,Me] = ke_solid10(Node,Section,Material,Options) computes element 
%   stiffness and mass matrix in the global coordinate system for a solid10 
%   element.
%
%   Node       Node definitions           [x y z] (10 * 3)
%              Nodes should have the following order:
%
%                    4
%                  / | \
%                 8  | 10
%                /   9   \
%               1--- |-7--3
%                \   |   /
%                 5  |  6
%                  \ | /
%                    2
%
%   Section    Section definition         []
%   Material   Material definition        [E nu rho]
%   Options    Struct containing optional parameters. Fields:
%      .nXi Number of Gauss integration points
%            {1 | 4 (default) | 5}
%   Ke         Element stiffness matrix (30 * 30)
%   Me         Element mass matrix (30 * 30)
%
%   See also KE_BEAM, ASMKM, KE_TRUSS.

if nargin<4, Options=[]; end
if ~isfield(Options,'nXi'), Options.nXi=4; end

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

if (nargout==2), rho=Material(3); end

% INITIALIZATION
Ke=zeros(30,30);
if nargout>1, Me = []; end
if (nargout==2),  Me=zeros(30,30); end

C=cmat_isotropic('3d',Section,Material);

% --- Integration points and weights---
nXi=Options.nXi;
[xi,H] = gaussqtet(nXi);

% --- Shape function ---
N=zeros(nXi,10);
N(:, 1)=(xi(:,1) + xi(:,2) + xi(:,3) - 1).*(2*xi(:,1) + 2*xi(:,2) + 2*xi(:,3) - 1); % N1
N(:, 2)=xi(:,1).*(2*xi(:,1) - 1);                                                   % N2
N(:, 3)=xi(:,2).*(2*xi(:,2) - 1);                                                   % N3
N(:, 4)=xi(:,3).*(2*xi(:,3) - 1);                                                   % N4
N(:, 5)=-xi(:,1).*(4*xi(:,1) + 4*xi(:,2) + 4*xi(:,3) - 4);                          % N5
N(:, 6)=4*xi(:,1).*xi(:,2);                                                         % N6
N(:, 7)=-4*xi(:,2).*(xi(:,1) + xi(:,2) + xi(:,3) - 1);                              % N7
N(:, 8)=-xi(:,3).*(4*xi(:,1) + 4*xi(:,2) + 4*xi(:,3) - 4);                          % N8
N(:, 9)=4*xi(:,1).*xi(:,3);                                                         % N9
N(:,10)=4*xi(:,2).*xi(:,3);                                                         % N10

% --- Shape function derivative ---
dN=zeros(nXi,30);
dN(:, 1)= 4*xi(:,1)+4*xi(:,2)+4*xi(:,3)-3;   % dN1/dxi1
dN(:, 2)= 4*xi(:,1)-1;                       % dN2/dxi1
dN(:, 3)= 0;                                 % dN3/dxi1
dN(:, 4)= 0;                                 % dN4/dxi1
dN(:, 5)= 4-4*xi(:,2)-4*xi(:,3)-8*xi(:,1);   % dN5/dxi1
dN(:, 6)= 4*xi(:,2);                         % dN6/dxi1
dN(:, 7)=-4*xi(:,2);                         % dN7/dxi1
dN(:, 8)=-4*xi(:,3);                         % dN8/dxi1
dN(:, 9)= 4*xi(:,3);                         % dN9/dxi1
dN(:,10)= 0;                                 % dN10/dxi1

dN(:,11)= 4*xi(:,1)+4*xi(:,2)+4*xi(:,3)-3;   % dN1/dxi2
dN(:,12)= 0;                                 % dN2/dxi2
dN(:,13)= 4*xi(:,2)-1;                       % dN3/dxi2
dN(:,14)= 0;                                 % dN4/dxi2
dN(:,15)=-4*xi(:,1);                         % dN5/dxi2
dN(:,16)= 4*xi(:,1);                         % dN6/dxi2
dN(:,17)= 4-8*xi(:,2)-4*xi(:,3)-4*xi(:,1);   % dN7/dxi2
dN(:,18)=-4*xi(:,3);                         % dN8/dxi2
dN(:,19)= 0;                                 % dN9/dxi2
dN(:,20)= 4*xi(:,3);                         % dN10/dxi2

dN(:,21)=4*xi(:,1)+4*xi(:,2)+4*xi(:,3)-3;    % dN1/dxi3
dN(:,22)=0;                                  % dN2/dxi3
dN(:,23)=0;                                  % dN3/dxi3
dN(:,24)=4*xi(:,3)-1;                        % dN4/dxi3
dN(:,25)=-4*xi(:,1);                         % dN5/dxi3
dN(:,26)=0;                                  % dN6/dxi3
dN(:,27)=-4*xi(:,2);                         % dN7/dxi3
dN(:,28)=4-4*xi(:,2)-8*xi(:,3)-4*xi(:,1);    % dN8/dxi3
dN(:,29)=4*xi(:,1);                          % dN9/dxi3
dN(:,30)=4*xi(:,2);                          % dN10/dxi3

% --- Jacobian ---
% Jac =  [dx/dxi1 dx/dxi2 dx/dxi3
%         dy/dxi1 dy/dxi2 dy/dxi3
%         dz/dxi1 dz/dxi2 dz/dxi3]
Jac=zeros(nXi,3,3);
for iXi=1:nXi
	for iNod=1:10
	  Jac(iXi,1,1)=Jac(iXi,1,1)+dN(iXi,iNod)*Node(iNod,1);  
	  Jac(iXi,1,2)=Jac(iXi,1,2)+dN(iXi,iNod)*Node(iNod,2);  
    Jac(iXi,1,3)=Jac(iXi,1,3)+dN(iXi,iNod)*Node(iNod,3);  
	  Jac(iXi,2,1)=Jac(iXi,2,1)+dN(iXi,10+iNod)*Node(iNod,1);
	  Jac(iXi,2,2)=Jac(iXi,2,2)+dN(iXi,10+iNod)*Node(iNod,2);
    Jac(iXi,2,3)=Jac(iXi,2,3)+dN(iXi,10+iNod)*Node(iNod,3);
    Jac(iXi,3,1)=Jac(iXi,3,1)+dN(iXi,20+iNod)*Node(iNod,1);
	  Jac(iXi,3,2)=Jac(iXi,3,2)+dN(iXi,20+iNod)*Node(iNod,2);
    Jac(iXi,3,3)=Jac(iXi,3,3)+dN(iXi,20+iNod)*Node(iNod,3);
  end
end

% --- Integration: loop over integration points ---
for iXi=1:nXi
  JacUtil=[Jac(iXi,1,1),Jac(iXi,1,2),Jac(iXi,1,3);Jac(iXi,2,1),Jac(iXi,2,2),Jac(iXi,2,3);Jac(iXi,3,1),Jac(iXi,3,2),Jac(iXi,3,3)];
  detJac=det(JacUtil);
  DN1= JacUtil\[dN(iXi, 1);dN(iXi,11);dN(iXi,21)];
  DN2= JacUtil\[dN(iXi, 2);dN(iXi,12);dN(iXi,22)];
  DN3= JacUtil\[dN(iXi, 3);dN(iXi,13);dN(iXi,23)];
  DN4= JacUtil\[dN(iXi, 4);dN(iXi,14);dN(iXi,24)];
  DN5= JacUtil\[dN(iXi, 5);dN(iXi,15);dN(iXi,25)];
  DN6= JacUtil\[dN(iXi, 6);dN(iXi,16);dN(iXi,26)];
  DN7= JacUtil\[dN(iXi, 7);dN(iXi,17);dN(iXi,27)];
  DN8= JacUtil\[dN(iXi, 8);dN(iXi,18);dN(iXi,28)];
  DN9= JacUtil\[dN(iXi, 9);dN(iXi,19);dN(iXi,29)];
  DN10=JacUtil\[dN(iXi,10);dN(iXi,20);dN(iXi,30)];
  
B=[DN1(1)  0        0      DN2(1)  0        0      DN3(1)  0        0      DN4(1)  0        0      DN5(1)  0        0      DN6(1)  0        0      DN7(1)  0        0      DN8(1)  0        0       DN9(1)  0        0      DN10(1)  0         0       
   0       DN1(2)   0      0       DN2(2)   0      0       DN3(2)   0      0       DN4(2)   0      0       DN5(2)   0      0       DN6(2)   0      0       DN7(2)   0      0       DN8(2)   0       0       DN9(2)   0      0        DN10(2)   0       
   0       0        DN1(3) 0       0        DN2(3) 0       0        DN3(3) 0       0        DN4(3) 0       0        DN5(3) 0       0        DN6(3) 0       0        DN7(3) 0       0        DN8(3)  0       0        DN9(3) 0        0         DN10(3) 
   DN1(2)  DN1(1)   0      DN2(2)  DN2(1)   0      DN3(2)  DN3(1)   0      DN4(2)  DN4(1)   0      DN5(2)  DN5(1)   0      DN6(2)  DN6(1)   0      DN7(2)  DN7(1)   0      DN8(2)  DN8(1)   0       DN9(2)  DN9(1)   0      DN10(2)  DN10(1)   0       
   0       DN1(3)   DN1(2) 0       DN2(3)   DN2(2) 0       DN3(3)   DN3(2) 0       DN4(3)   DN4(2) 0       DN5(3)   DN5(2) 0       DN6(3)   DN6(2) 0       DN7(3)   DN7(2) 0       DN8(3)   DN8(2)  0       DN9(3)   DN9(2) 0        DN10(3)   DN10(2) 
   DN1(3)  0        DN1(1) DN2(3)  0        DN2(1) DN3(3)  0        DN3(1) DN4(3)  0        DN4(1) DN5(3)  0        DN5(1) DN6(3)  0        DN6(1) DN7(3)  0        DN7(1) DN8(3)  0        DN8(1)  DN9(3)  0        DN9(1) DN10(3)  0         DN10(1) ];

  Ke=Ke+H(iXi)*detJac*(B.'*C*B);
  if (nargout==2)
    Nxi=[N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        0        N(iXi,9) 0        0        N(iXi,10) 0         0        
         0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        0        N(iXi,9) 0        0         N(iXi,10) 0        
         0        0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        0        N(iXi,9) 0         0         N(iXi,10)];
    Me=Me+H(iXi)*detJac*rho*(Nxi.'*Nxi);
  end
end

