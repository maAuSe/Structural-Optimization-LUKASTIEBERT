function [Ke,Me,dKedx]=ke_solid4(Node,Section,Material,Options,dNodedx,dSectiondx)
%KE_SOLID4   Compute the element stiffness and mass matrix for a solid4 element.
%
%   [Ke,Me] = ke_solid4(Node,Section,Material,Options) computes element 
%   stiffness and mass matrix in the global coordinate system for a solid4 
%   element.
%
%   Node       Node definitions           [x y z] (4 * 3)
%              Nodes should have the following order:
%
%                    4
%                  / | \
%                 /  |  \
%                /   |   \
%               1--- |----3
%                \   |   /
%                 \  |  /
%                  \ | /
%                    2
%
%   Section    Section definition         []
%   Material   Material definition        [E nu rho]
%   Options    Element options struct. Fields: []
%   Ke         Element stiffness matrix (12 * 12)
%   Me         Element mass matrix (12 * 12)
%
%   See also KE_BEAM, ASMKM, KE_TRUSS.

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
if (nargout==2), rho=Material(3); end

[~,C_lambda,C_mu]=cmat_isotropic('3d',Section,Material);

% INITIALIZATION
Ke=zeros(12,12);
if nargout>1, Me = []; end
if (nargout==2),  Me=zeros(12,12); end

%===============================================================================
% SELECTIVE INTEGRATION: INTEGRATION OF MU-TERM
%===============================================================================

% --- Integration points and weights---
nXi=4;
[xi,H] = gaussqtet(nXi);

% --- Shape function ---
N=zeros(nXi,4);
N(:,1)=1-xi(:,1)- xi(:,2)-xi(:,3); % N1(xi1,xi2,xi3)  
N(:,2)=xi(:,1);                    % N2(xi1,xi2,xi3)  
N(:,3)=xi(:,2);                    % N3(xi1,xi2,xi3)  
N(:,4)=xi(:,3);                    % N4(xi1,xi2,xi3)  

% --- Shape function derivative ---
dN=zeros(nXi,12);
dN(:, 1)=-1;                      % dN1(xi1,xi2,xi3)/dxi1
dN(:, 2)= 1;                      % dN2(xi1,xi2,xi3)/dxi1
dN(:, 3)= 0;                      % dN3(xi1,xi2,xi3)/dxi1
dN(:, 4)= 0;                      % dN4(xi1,xi2,xi3)/dxi1
dN(:, 5)=-1;                      % dN1(xi1,xi2,xi3)/dxi2
dN(:, 6)= 0;                      % dN2(xi1,xi2,xi3)/dxi2
dN(:, 7)= 1;                      % dN3(xi1,xi2,xi3)/dxi2
dN(:, 8)= 0;                      % dN4(xi1,xi2,xi3)/dxi2
dN(:, 9)=-1;                      % dN1(xi1,xi2,xi3)/dxi3
dN(:,10)= 0;                      % dN2(xi1,xi2,xi3)/dxi3
dN(:,11)= 0;                      % dN3(xi1,xi2,xi3)/dxi3
dN(:,12)= 1;                      % dN4(xi1,xi2,xi3)/dxi3

% --- Jacobian ---
% Jac =  [dx/dxi1 dx/dxi2 dx/dxi3
%         dy/dxi1 dy/dxi2 dy/dxi3
%         dz/dxi1 dz/dxi2 dz/dxi3]
Jac=zeros(nXi,3,3);
for iXi=1:nXi
	for iNod=1:4
	  Jac(iXi,1,1)=Jac(iXi,1,1)+dN(iXi,iNod)   * Node(iNod,1);  
	  Jac(iXi,1,2)=Jac(iXi,1,2)+dN(iXi,iNod)   * Node(iNod,2);  
    Jac(iXi,1,3)=Jac(iXi,1,3)+dN(iXi,iNod)   * Node(iNod,3);  
	  Jac(iXi,2,1)=Jac(iXi,2,1)+dN(iXi,4+iNod) * Node(iNod,1);
	  Jac(iXi,2,2)=Jac(iXi,2,2)+dN(iXi,4+iNod) * Node(iNod,2);
    Jac(iXi,2,3)=Jac(iXi,2,3)+dN(iXi,4+iNod) * Node(iNod,3);
    Jac(iXi,3,1)=Jac(iXi,3,1)+dN(iXi,8+iNod) * Node(iNod,1);
	  Jac(iXi,3,2)=Jac(iXi,3,2)+dN(iXi,8+iNod) * Node(iNod,2);
    Jac(iXi,3,3)=Jac(iXi,3,3)+dN(iXi,8+iNod) * Node(iNod,3);
  end
end

% --- Integration: loop over integration points ---
for iXi=1:nXi
  JacUtil=[Jac(iXi,1,1),Jac(iXi,1,2),Jac(iXi,1,3);Jac(iXi,2,1),Jac(iXi,2,2),Jac(iXi,2,3);Jac(iXi,3,1),Jac(iXi,3,2),Jac(iXi,3,3)];
  detJac=det(JacUtil);
  DN1= JacUtil\[dN(iXi,1);dN(iXi,5);dN(iXi,9)];     % [dN1/dx; dN1/dy; dN1/dz]
  DN2= JacUtil\[dN(iXi,2);dN(iXi,6);dN(iXi,10)];    % [dN2/dx; dN2/dy; dN2/dz]
  DN3= JacUtil\[dN(iXi,3);dN(iXi,7);dN(iXi,11)];    % [dN3/dx; dN3/dy; dN3/dz]
  DN4= JacUtil\[dN(iXi,4);dN(iXi,8);dN(iXi,12)];    % [dN4/dx; dN4/dy; dN4/dz]

  B=[DN1(1)  0        0      DN2(1)  0        0      DN3(1)  0        0      DN4(1)  0        0        
     0       DN1(2)   0      0       DN2(2)   0      0       DN3(2)   0      0       DN4(2)   0        
     0       0        DN1(3) 0       0        DN2(3) 0       0        DN3(3) 0       0        DN4(3)   
     DN1(2)  DN1(1)   0      DN2(2)  DN2(1)   0      DN3(2)  DN3(1)   0      DN4(2)  DN4(1)   0        
     0       DN1(3)   DN1(2) 0       DN2(3)   DN2(2) 0       DN3(3)   DN3(2) 0       DN4(3)   DN4(2)   
     DN1(3)  0        DN1(1) DN2(3)  0        DN2(1) DN3(3)  0        DN3(1) DN4(3)  0        DN4(1)];

  Ke=Ke+H(iXi)*detJac*(B.'*C_mu*B);
  if (nargout==2)
      Nxi=[N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        
           0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        
           0        0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4)];
      Me=Me+H(iXi)*detJac*rho*(Nxi.'*Nxi);
  end
end


%===============================================================================
% SELECTIVE INTEGRATION: INTEGRATION OF LAMBDA-TERM
%===============================================================================

% --- Integration points and weights---
nXi=1;
[xi,H] = gaussqtet(nXi);

% --- Shape function ---
N=zeros(nXi,4);
N(:,1)=1-xi(:,1)- xi(:,2)-xi(:,3); % N1(xi1,xi2,xi3)  
N(:,2)=xi(:,1);                    % N2(xi1,xi2,xi3)  
N(:,3)=xi(:,2);                    % N3(xi1,xi2,xi3)  
N(:,4)=xi(:,3);                    % N4(xi1,xi2,xi3)  

% --- Shape function derivative ---
dN=zeros(nXi,12);
dN(:, 1)=-1;                      % dN1(xi1,xi2,xi3)/dxi1
dN(:, 2)= 1;                      % dN2(xi1,xi2,xi3)/dxi1
dN(:, 3)= 0;                      % dN3(xi1,xi2,xi3)/dxi1
dN(:, 4)= 0;                      % dN4(xi1,xi2,xi3)/dxi1
dN(:, 5)=-1;                      % dN1(xi1,xi2,xi3)/dxi2
dN(:, 6)= 0;                      % dN2(xi1,xi2,xi3)/dxi2
dN(:, 7)= 1;                      % dN3(xi1,xi2,xi3)/dxi2
dN(:, 8)= 0;                      % dN4(xi1,xi2,xi3)/dxi2
dN(:, 9)=-1;                      % dN1(xi1,xi2,xi3)/dxi3
dN(:,10)= 0;                      % dN2(xi1,xi2,xi3)/dxi3
dN(:,11)= 0;                      % dN3(xi1,xi2,xi3)/dxi3
dN(:,12)= 1;                      % dN4(xi1,xi2,xi3)/dxi3

% --- Jacobian ---
% Jac =  [dx/dxi1 dx/dxi2 dx/dxi3
%         dy/dxi1 dy/dxi2 dy/dxi3
%         dz/dxi1 dz/dxi2 dz/dxi3]
Jac=zeros(nXi,3,3);
for iXi=1:nXi
	for iNod=1:4
	  Jac(iXi,1,1)=Jac(iXi,1,1)+dN(iXi,iNod)   * Node(iNod,1);  
	  Jac(iXi,1,2)=Jac(iXi,1,2)+dN(iXi,iNod)   * Node(iNod,2);  
    Jac(iXi,1,3)=Jac(iXi,1,3)+dN(iXi,iNod)   * Node(iNod,3);  
	  Jac(iXi,2,1)=Jac(iXi,2,1)+dN(iXi,4+iNod) * Node(iNod,1);
	  Jac(iXi,2,2)=Jac(iXi,2,2)+dN(iXi,4+iNod) * Node(iNod,2);
    Jac(iXi,2,3)=Jac(iXi,2,3)+dN(iXi,4+iNod) * Node(iNod,3);
    Jac(iXi,3,1)=Jac(iXi,3,1)+dN(iXi,8+iNod) * Node(iNod,1);
	  Jac(iXi,3,2)=Jac(iXi,3,2)+dN(iXi,8+iNod) * Node(iNod,2);
    Jac(iXi,3,3)=Jac(iXi,3,3)+dN(iXi,8+iNod) * Node(iNod,3);
  end
end

% --- Integration: loop over integration points ---
for iXi=1:nXi
  JacUtil=[Jac(iXi,1,1),Jac(iXi,1,2),Jac(iXi,1,3);Jac(iXi,2,1),Jac(iXi,2,2),Jac(iXi,2,3);Jac(iXi,3,1),Jac(iXi,3,2),Jac(iXi,3,3)];
  detJac=det(JacUtil);
  DN1= JacUtil\[dN(iXi,1);dN(iXi,5);dN(iXi,9)];     % [dN1/dx; dN1/dy; dN1/dz]
  DN2= JacUtil\[dN(iXi,2);dN(iXi,6);dN(iXi,10)];    % [dN2/dx; dN2/dy; dN2/dz]
  DN3= JacUtil\[dN(iXi,3);dN(iXi,7);dN(iXi,11)];    % [dN3/dx; dN3/dy; dN3/dz]
  DN4= JacUtil\[dN(iXi,4);dN(iXi,8);dN(iXi,12)];    % [dN4/dx; dN4/dy; dN4/dz]

  B=[DN1(1)  0        0      DN2(1)  0        0      DN3(1)  0        0      DN4(1)  0        0        
     0       DN1(2)   0      0       DN2(2)   0      0       DN3(2)   0      0       DN4(2)   0        
     0       0        DN1(3) 0       0        DN2(3) 0       0        DN3(3) 0       0        DN4(3)   
     DN1(2)  DN1(1)   0      DN2(2)  DN2(1)   0      DN3(2)  DN3(1)   0      DN4(2)  DN4(1)   0        
     0       DN1(3)   DN1(2) 0       DN2(3)   DN2(2) 0       DN3(3)   DN3(2) 0       DN4(3)   DN4(2)   
     DN1(3)  0        DN1(1) DN2(3)  0        DN2(1) DN3(3)  0        DN3(1) DN4(3)  0        DN4(1)];

  Ke=Ke+H(iXi)*detJac*(B.'*C_lambda*B);
end