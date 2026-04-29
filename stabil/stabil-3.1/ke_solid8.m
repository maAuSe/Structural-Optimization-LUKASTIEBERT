function [Ke,Me,dKedx]=ke_solid8(Node,Section,Material,Options,dNodedx,dSectiondx)
%KE_SOLID8   Compute the element stiffness and mass matrix for a solid8 element.
%
%   [Ke,Me] = ke_solid8(Node,Section,Material,Options)
%   [Ke,Me] = ke_solid8(Node,Section,Material,Options)
%   computes element stiffness and mass matrix in the global coordinate system 
%                                                             for a solid8 element.
%
%   Node       Node definitions           [x y z] (8 * 3)
%              Nodes should have the following order:
%                  8---------7
%                 /|        /|
%                / |       / |
%               /  |      /  |
%              5---------6   |
%              |   4-----|---3
%              |  /      |  /
%              | /       | /
%              |/        |/
%              1---------2
%
%   Section    Section definition         []
%   Material   Material definition        [E nu rho]
%   Options    Element options struct. Fields: []
%   Ke         Element stiffness matrix (24 * 24)
%   Me         Element mass matrix (24 * 24)
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

% MATERIAL PROPERTIES
if (nargout==2), rho=Material(3); end

[~,C_lambda,C_mu]=cmat_isotropic('3d',Section,Material);

% INITIALIZATION
Ke=zeros(24,24);
if nargout>1, Me = []; end
if (nargout==2),  Me=zeros(24,24); end

%===============================================================================
% SELECTIVE INTEGRATION: INTEGRATION OF MU-TERM
%===============================================================================

% --- Integration points and weights---
nXi1D=2;
xi1D=[5.773502691896258e-001;-5.773502691896258e-001];
H1D =[1;1];

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
N=zeros(nXi,8);
N(:,1)=1/8.*(1-xi(:,1)).*(1-xi(:,2)).*(1-xi(:,3));  % N1(xi1,xi2,xi3)  
N(:,2)=1/8.*(1+xi(:,1)).*(1-xi(:,2)).*(1-xi(:,3));  % N2(xi1,xi2,xi3)  
N(:,3)=1/8.*(1+xi(:,1)).*(1+xi(:,2)).*(1-xi(:,3));  % N3(xi1,xi2,xi3)  
N(:,4)=1/8.*(1-xi(:,1)).*(1+xi(:,2)).*(1-xi(:,3));  % N4(xi1,xi2,xi3)  
N(:,5)=1/8.*(1-xi(:,1)).*(1-xi(:,2)).*(1+xi(:,3));  % N5(xi1,xi2,xi3)  
N(:,6)=1/8.*(1+xi(:,1)).*(1-xi(:,2)).*(1+xi(:,3));  % N6(xi1,xi2,xi3)  
N(:,7)=1/8.*(1+xi(:,1)).*(1+xi(:,2)).*(1+xi(:,3));  % N7(xi1,xi2,xi3)  
N(:,8)=1/8.*(1-xi(:,1)).*(1+xi(:,2)).*(1+xi(:,3));  % N8(xi1,xi2,xi3)  

% --- Shape function derivative ---
dN=zeros(nXi,24);
dN(:,1) =-(xi(:,2)-1).*(xi(:,3)-1)/8;        % dN1(xi1,xi2,xi3)/dxi1
dN(:,2) = (xi(:,2)-1).*(xi(:,3)-1)/8;        % dN2(xi1,xi2,xi3)/dxi1
dN(:,3) =-(xi(:,2)+1).*(xi(:,3)-1)/8;        % dN3(xi1,xi2,xi3)/dxi1
dN(:,4) = (xi(:,2)+1).*(xi(:,3)-1)/8;        % dN4(xi1,xi2,xi3)/dxi1
dN(:,5) = (xi(:,2)-1).*(xi(:,3)+1)/8;        % dN5(xi1,xi2,xi3)/dxi1
dN(:,6) =-(xi(:,2)-1).*(xi(:,3)+1)/8;        % dN6(xi1,xi2,xi3)/dxi1
dN(:,7) = (xi(:,2)+1).*(xi(:,3)+1)/8;        % dN7(xi1,xi2,xi3)/dxi1
dN(:,8) =-(xi(:,2)+1).*(xi(:,3)+1)/8;        % dN8(xi1,xi2,xi3)/dxi1
dN(:,9) =-(xi(:,1)-1).*(xi(:,3)-1)/8;        % dN1(xi1,xi2,xi3)/dxi2
dN(:,10)= (xi(:,1)+1).*(xi(:,3)-1)/8;        % dN2(xi1,xi2,xi3)/dxi2
dN(:,11)=-(xi(:,1)+1).*(xi(:,3)-1)/8;        % dN3(xi1,xi2,xi3)/dxi2
dN(:,12)= (xi(:,1)-1).*(xi(:,3)-1)/8;        % dN4(xi1,xi2,xi3)/dxi2
dN(:,13)= (xi(:,1)-1).*(xi(:,3)+1)/8;        % dN5(xi1,xi2,xi3)/dxi2
dN(:,14)=-(xi(:,1)+1).*(xi(:,3)+1)/8;        % dN6(xi1,xi2,xi3)/dxi2
dN(:,15)= (xi(:,1)+1).*(xi(:,3)+1)/8;        % dN7(xi1,xi2,xi3)/dxi2
dN(:,16)=-(xi(:,1)-1).*(xi(:,3)+1)/8;        % dN8(xi1,xi2,xi3)/dxi2
dN(:,17)=-(xi(:,1)-1).*(xi(:,2)-1)/8;        % dN1(xi1,xi2,xi3)/dxi3
dN(:,18)= (xi(:,1)+1).*(xi(:,2)-1)/8;        % dN2(xi1,xi2,xi3)/dxi3
dN(:,19)=-(xi(:,1)+1).*(xi(:,2)+1)/8;        % dN3(xi1,xi2,xi3)/dxi3
dN(:,20)= (xi(:,1)-1).*(xi(:,2)+1)/8;        % dN4(xi1,xi2,xi3)/dxi3
dN(:,21)= (xi(:,1)-1).*(xi(:,2)-1)/8;        % dN5(xi1,xi2,xi3)/dxi3
dN(:,22)=-(xi(:,1)+1).*(xi(:,2)-1)/8;        % dN6(xi1,xi2,xi3)/dxi3
dN(:,23)= (xi(:,1)+1).*(xi(:,2)+1)/8;        % dN7(xi1,xi2,xi3)/dxi3
dN(:,24)=-(xi(:,1)-1).*(xi(:,2)+1)/8;        % dN8(xi1,xi2,xi3)/dxi3

% --- Jacobian ---
% Jac =  [dx/dxi1 dx/dxi2 dx/dxi3
%         dy/dxi1 dy/dxi2 dy/dxi3
%         dz/dxi1 dz/dxi2 dz/dxi3]
Jac=zeros(nXi,3,3);
for iXi=1:nXi
	for iNod=1:8
	  Jac(iXi,1,1)=Jac(iXi,1,1)+dN(iXi,iNod)    * Node(iNod,1);  
	  Jac(iXi,1,2)=Jac(iXi,1,2)+dN(iXi,iNod)    * Node(iNod,2);  
    Jac(iXi,1,3)=Jac(iXi,1,3)+dN(iXi,iNod)    * Node(iNod,3);  
	  Jac(iXi,2,1)=Jac(iXi,2,1)+dN(iXi,8+iNod)  * Node(iNod,1);
	  Jac(iXi,2,2)=Jac(iXi,2,2)+dN(iXi,8+iNod)  * Node(iNod,2);
    Jac(iXi,2,3)=Jac(iXi,2,3)+dN(iXi,8+iNod)  * Node(iNod,3);
    Jac(iXi,3,1)=Jac(iXi,3,1)+dN(iXi,16+iNod) * Node(iNod,1);
	  Jac(iXi,3,2)=Jac(iXi,3,2)+dN(iXi,16+iNod) * Node(iNod,2);
    Jac(iXi,3,3)=Jac(iXi,3,3)+dN(iXi,16+iNod) * Node(iNod,3);
  end
end

% --- Integration: loop over integration points ---
for iXi=1:nXi
  JacUtil=[Jac(iXi,1,1),Jac(iXi,1,2),Jac(iXi,1,3);Jac(iXi,2,1),Jac(iXi,2,2),Jac(iXi,2,3);Jac(iXi,3,1),Jac(iXi,3,2),Jac(iXi,3,3)];
  detJac=det(JacUtil);
  DN1= JacUtil\[dN(iXi,1);dN(iXi,9) ;dN(iXi,17)];    % [dN1/dx; dN1/dy; dN1/dz]
  DN2= JacUtil\[dN(iXi,2);dN(iXi,10);dN(iXi,18)];    % [dN2/dx; dN2/dy; dN2/dz]
  DN3= JacUtil\[dN(iXi,3);dN(iXi,11);dN(iXi,19)];    % [dN3/dx; dN3/dy; dN3/dz]
  DN4= JacUtil\[dN(iXi,4);dN(iXi,12);dN(iXi,20)];    % [dN4/dx; dN4/dy; dN4/dz]
  DN5= JacUtil\[dN(iXi,5);dN(iXi,13);dN(iXi,21)];    % [dN5/dx; dN5/dy; dN5/dz]
  DN6= JacUtil\[dN(iXi,6);dN(iXi,14);dN(iXi,22)];    % [dN6/dx; dN6/dy; dN6/dz]
  DN7= JacUtil\[dN(iXi,7);dN(iXi,15);dN(iXi,23)];    % [dN7/dx; dN7/dy; dN7/dz]
  DN8= JacUtil\[dN(iXi,8);dN(iXi,16);dN(iXi,24)];    % [dN8/dx; dN8/dy; dN8/dz]

  B=[DN1(1)  0        0      DN2(1)  0        0      DN3(1)  0        0      DN4(1)  0        0      DN5(1)  0        0      DN6(1)  0        0      DN7(1)  0        0      DN8(1)  0        0        
     0       DN1(2)   0      0       DN2(2)   0      0       DN3(2)   0      0       DN4(2)   0      0       DN5(2)   0      0       DN6(2)   0      0       DN7(2)   0      0       DN8(2)   0        
     0       0        DN1(3) 0       0        DN2(3) 0       0        DN3(3) 0       0        DN4(3) 0       0        DN5(3) 0       0        DN6(3) 0       0        DN7(3) 0       0        DN8(3)   
     DN1(2)  DN1(1)   0      DN2(2)  DN2(1)   0      DN3(2)  DN3(1)   0      DN4(2)  DN4(1)   0      DN5(2)  DN5(1)   0      DN6(2)  DN6(1)   0      DN7(2)  DN7(1)   0      DN8(2)  DN8(1)   0        
     0       DN1(3)   DN1(2) 0       DN2(3)   DN2(2) 0       DN3(3)   DN3(2) 0       DN4(3)   DN4(2) 0       DN5(3)   DN5(2) 0       DN6(3)   DN6(2) 0       DN7(3)   DN7(2) 0       DN8(3)   DN8(2)   
     DN1(3)  0        DN1(1) DN2(3)  0        DN2(1) DN3(3)  0        DN3(1) DN4(3)  0        DN4(1) DN5(3)  0        DN5(1) DN6(3)  0        DN6(1) DN7(3)  0        DN7(1) DN8(3)  0        DN8(1)];

  Nxi=[N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        0        
       0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        
       0        0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8)];
  Ke=Ke+H(iXi)*detJac*(B.'*C_mu*B);
  if (nargout==2)
    Me=Me+H(iXi)*detJac*rho*(Nxi.'*Nxi);
  end
end


%===============================================================================
% SELECTIVE INTEGRATION: INTEGRATION OF LAMBDA-TERM
%===============================================================================

nXi1D=1; % single integration point used
xi1D=0;
H1D =2;

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
N=zeros(nXi,8);
N(:,1)=1/8.*(1-xi(:,1)).*(1-xi(:,2)).*(1-xi(:,3));  % N1(xi1,xi2,xi3)  
N(:,2)=1/8.*(1+xi(:,1)).*(1-xi(:,2)).*(1-xi(:,3));  % N2(xi1,xi2,xi3)  
N(:,3)=1/8.*(1+xi(:,1)).*(1+xi(:,2)).*(1-xi(:,3));  % N3(xi1,xi2,xi3)  
N(:,4)=1/8.*(1-xi(:,1)).*(1+xi(:,2)).*(1-xi(:,3));  % N4(xi1,xi2,xi3)  
N(:,5)=1/8.*(1-xi(:,1)).*(1-xi(:,2)).*(1+xi(:,3));  % N5(xi1,xi2,xi3)  
N(:,6)=1/8.*(1+xi(:,1)).*(1-xi(:,2)).*(1+xi(:,3));  % N6(xi1,xi2,xi3)  
N(:,7)=1/8.*(1+xi(:,1)).*(1+xi(:,2)).*(1+xi(:,3));  % N7(xi1,xi2,xi3)  
N(:,8)=1/8.*(1-xi(:,1)).*(1+xi(:,2)).*(1+xi(:,3));  % N8(xi1,xi2,xi3)  

% --- Shape function derivative ---
dN=zeros(nXi,24);
dN(:,1) =-(xi(:,2)-1).*(xi(:,3)-1)/8;        % dN1(xi1,xi2,xi3)/dxi1
dN(:,2) = (xi(:,2)-1).*(xi(:,3)-1)/8;        % dN2(xi1,xi2,xi3)/dxi1
dN(:,3) =-(xi(:,2)+1).*(xi(:,3)-1)/8;        % dN3(xi1,xi2,xi3)/dxi1
dN(:,4) = (xi(:,2)+1).*(xi(:,3)-1)/8;        % dN4(xi1,xi2,xi3)/dxi1
dN(:,5) = (xi(:,2)-1).*(xi(:,3)+1)/8;        % dN5(xi1,xi2,xi3)/dxi1
dN(:,6) =-(xi(:,2)-1).*(xi(:,3)+1)/8;        % dN6(xi1,xi2,xi3)/dxi1
dN(:,7) = (xi(:,2)+1).*(xi(:,3)+1)/8;        % dN7(xi1,xi2,xi3)/dxi1
dN(:,8) =-(xi(:,2)+1).*(xi(:,3)+1)/8;        % dN8(xi1,xi2,xi3)/dxi1
dN(:,9) =-(xi(:,1)-1).*(xi(:,3)-1)/8;        % dN1(xi1,xi2,xi3)/dxi2
dN(:,10)= (xi(:,1)+1).*(xi(:,3)-1)/8;        % dN2(xi1,xi2,xi3)/dxi2
dN(:,11)=-(xi(:,1)+1).*(xi(:,3)-1)/8;        % dN3(xi1,xi2,xi3)/dxi2
dN(:,12)= (xi(:,1)-1).*(xi(:,3)-1)/8;        % dN4(xi1,xi2,xi3)/dxi2
dN(:,13)= (xi(:,1)-1).*(xi(:,3)+1)/8;        % dN5(xi1,xi2,xi3)/dxi2
dN(:,14)=-(xi(:,1)+1).*(xi(:,3)+1)/8;        % dN6(xi1,xi2,xi3)/dxi2
dN(:,15)= (xi(:,1)+1).*(xi(:,3)+1)/8;        % dN7(xi1,xi2,xi3)/dxi2
dN(:,16)=-(xi(:,1)-1).*(xi(:,3)+1)/8;        % dN8(xi1,xi2,xi3)/dxi2
dN(:,17)=-(xi(:,1)-1).*(xi(:,2)-1)/8;        % dN1(xi1,xi2,xi3)/dxi3
dN(:,18)= (xi(:,1)+1).*(xi(:,2)-1)/8;        % dN2(xi1,xi2,xi3)/dxi3
dN(:,19)=-(xi(:,1)+1).*(xi(:,2)+1)/8;        % dN3(xi1,xi2,xi3)/dxi3
dN(:,20)= (xi(:,1)-1).*(xi(:,2)+1)/8;        % dN4(xi1,xi2,xi3)/dxi3
dN(:,21)= (xi(:,1)-1).*(xi(:,2)-1)/8;        % dN5(xi1,xi2,xi3)/dxi3
dN(:,22)=-(xi(:,1)+1).*(xi(:,2)-1)/8;        % dN6(xi1,xi2,xi3)/dxi3
dN(:,23)= (xi(:,1)+1).*(xi(:,2)+1)/8;        % dN7(xi1,xi2,xi3)/dxi3
dN(:,24)=-(xi(:,1)-1).*(xi(:,2)+1)/8;        % dN8(xi1,xi2,xi3)/dxi3

% --- Jacobian ---
% Jac =  [dx/dxi1 dx/dxi2 dx/dxi3
%         dy/dxi1 dy/dxi2 dy/dxi3
%         dz/dxi1 dz/dxi2 dz/dxi3]
Jac=zeros(nXi,3,3);
for iXi=1:nXi
	for iNod=1:8
	  Jac(iXi,1,1)=Jac(iXi,1,1)+dN(iXi,iNod)    * Node(iNod,1);  
	  Jac(iXi,1,2)=Jac(iXi,1,2)+dN(iXi,iNod)    * Node(iNod,2);  
    Jac(iXi,1,3)=Jac(iXi,1,3)+dN(iXi,iNod)    * Node(iNod,3);  
	  Jac(iXi,2,1)=Jac(iXi,2,1)+dN(iXi,8+iNod)  * Node(iNod,1);
	  Jac(iXi,2,2)=Jac(iXi,2,2)+dN(iXi,8+iNod)  * Node(iNod,2);
    Jac(iXi,2,3)=Jac(iXi,2,3)+dN(iXi,8+iNod)  * Node(iNod,3);
    Jac(iXi,3,1)=Jac(iXi,3,1)+dN(iXi,16+iNod) * Node(iNod,1);
	  Jac(iXi,3,2)=Jac(iXi,3,2)+dN(iXi,16+iNod) * Node(iNod,2);
    Jac(iXi,3,3)=Jac(iXi,3,3)+dN(iXi,16+iNod) * Node(iNod,3);
  end
end

% --- Integration: loop over integration points ---
for iXi=1:nXi
  JacUtil=[Jac(iXi,1,1),Jac(iXi,1,2),Jac(iXi,1,3);Jac(iXi,2,1),Jac(iXi,2,2),Jac(iXi,2,3);Jac(iXi,3,1),Jac(iXi,3,2),Jac(iXi,3,3)];
  detJac=det(JacUtil);
  DN1= JacUtil\[dN(iXi,1);dN(iXi,9) ;dN(iXi,17)];    % [dN1/dx; dN1/dy; dN1/dz]
  DN2= JacUtil\[dN(iXi,2);dN(iXi,10);dN(iXi,18)];    % [dN2/dx; dN2/dy; dN2/dz]
  DN3= JacUtil\[dN(iXi,3);dN(iXi,11);dN(iXi,19)];    % [dN3/dx; dN3/dy; dN3/dz]
  DN4= JacUtil\[dN(iXi,4);dN(iXi,12);dN(iXi,20)];    % [dN4/dx; dN4/dy; dN4/dz]
  DN5= JacUtil\[dN(iXi,5);dN(iXi,13);dN(iXi,21)];    % [dN5/dx; dN5/dy; dN5/dz]
  DN6= JacUtil\[dN(iXi,6);dN(iXi,14);dN(iXi,22)];    % [dN6/dx; dN6/dy; dN6/dz]
  DN7= JacUtil\[dN(iXi,7);dN(iXi,15);dN(iXi,23)];    % [dN7/dx; dN7/dy; dN7/dz]
  DN8= JacUtil\[dN(iXi,8);dN(iXi,16);dN(iXi,24)];    % [dN8/dx; dN8/dy; dN8/dz]

  B=[DN1(1)  0        0      DN2(1)  0        0      DN3(1)  0        0      DN4(1)  0        0      DN5(1)  0        0      DN6(1)  0        0      DN7(1)  0        0      DN8(1)  0        0        
     0       DN1(2)   0      0       DN2(2)   0      0       DN3(2)   0      0       DN4(2)   0      0       DN5(2)   0      0       DN6(2)   0      0       DN7(2)   0      0       DN8(2)   0        
     0       0        DN1(3) 0       0        DN2(3) 0       0        DN3(3) 0       0        DN4(3) 0       0        DN5(3) 0       0        DN6(3) 0       0        DN7(3) 0       0        DN8(3)   
     DN1(2)  DN1(1)   0      DN2(2)  DN2(1)   0      DN3(2)  DN3(1)   0      DN4(2)  DN4(1)   0      DN5(2)  DN5(1)   0      DN6(2)  DN6(1)   0      DN7(2)  DN7(1)   0      DN8(2)  DN8(1)   0        
     0       DN1(3)   DN1(2) 0       DN2(3)   DN2(2) 0       DN3(3)   DN3(2) 0       DN4(3)   DN4(2) 0       DN5(3)   DN5(2) 0       DN6(3)   DN6(2) 0       DN7(3)   DN7(2) 0       DN8(3)   DN8(2)   
     DN1(3)  0        DN1(1) DN2(3)  0        DN2(1) DN3(3)  0        DN3(1) DN4(3)  0        DN4(1) DN5(3)  0        DN5(1) DN6(3)  0        DN6(1) DN7(3)  0        DN7(1) DN8(3)  0        DN8(1)];

  Ke=Ke+H(iXi)*detJac*(B.'*C_lambda*B);
end