function [Ke,Me,dKedx]=ke_plane8(Node,Section,Material,Options,dNodedx,dSectiondx)

%KE_PLANE8   plane element stiffness and mass matrix in global coordinate system.
%
%   [Ke,Me] = ke_plane8(Node,Section,Material,Options) returns the element
%   stiffness and mass matrix in the global coordinate system for a 8-node
%   plane element. Plane8 only operates in the 2D xy-plane so that
%   z-coordinates should be equal to zero.
%
%   Node       Node definitions           [x y z] (8 * 3)
%   Section    Section definition         [h] (only used in plane stress)
%   Material   Material definition        [E nu rho]
%   Options    Struct containing optional parameters. Fields:
%      .problem Plane stress or plane strain
%               {'2dstress' (default) | '2dstrain' | 'axisym'}
%   Ke         Element stiffness matrix (16 * 16)
%   Me         Element mass matrix (16 * 16)

% Options
if nargin<4, Options=[]; end
if ~isfield(Options,'problem'), Options.problem='2dstress'; end

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

% Constitutive matrix
C=cmat_isotropic(Options.problem,Section,Material);

% --- Integration points and weights---
xi1D=[0.906179845938664;0.538469310105683;0;-0.538469310105683;-0.906179845938664];
H1D=[0.236926885056189;0.478628670499366;0.568888888888889;0.478628670499366;0.236926885056189];
nXi1D=numel(xi1D);
nXi=nXi1D^2;
xi=zeros(nXi,2);
H=zeros(nXi,1);
for iXi=1:nXi1D
    for jXi=1:nXi1D
        xi(nXi1D*(iXi-1)+jXi,1)= xi1D(iXi);             % xi
        xi(nXi1D*(iXi-1)+jXi,2)= xi1D(jXi);             % eta
        H(nXi1D*(iXi-1)+jXi)=H1D(iXi)*H1D(jXi);
    end
end

% --- Shape function ---
N=zeros(nXi,8);
N(:,1)=0.25*(1-xi(:,1)).*(1-xi(:,2)).*(-1-xi(:,1)-xi(:,2));   % N1(xi,eta)
N(:,2)=0.25*(1+xi(:,1)).*(1-xi(:,2)).*(-1+xi(:,1)-xi(:,2));   % N2(xi,eta)
N(:,3)=0.25*(1+xi(:,1)).*(1+xi(:,2)).*(-1+xi(:,1)+xi(:,2));   % N3(xi,eta)
N(:,4)=0.25*(1-xi(:,1)).*(1+xi(:,2)).*(-1-xi(:,1)+xi(:,2));   % N4(xi,eta)
N(:,5)=0.50*(1-xi(:,2)).*(1-xi(:,1).*xi(:,1));                % N5(xi,eta)
N(:,6)=0.50*(1+xi(:,1)).*(1-xi(:,2).*xi(:,2));                % N6(xi,eta)
N(:,7)=0.50*(1+xi(:,2)).*(1-xi(:,1).*xi(:,1));                % N7(xi,eta)
N(:,8)=0.50*(1-xi(:,1)).*(1-xi(:,2).*xi(:,2));                % N8(xi,eta)

% --- Shape function derivative ---
dN=zeros(nXi,16);
dN(:,1) = 0.25*(1.0-xi(:,2)).*(2.0*xi(:,1)+xi(:,2));         % dN1(xi,eta)/dxi
dN(:,2) = 0.25*(1.0-xi(:,2)).*(2.0*xi(:,1)-xi(:,2));         % dN2(xi,eta)/dxi
dN(:,3) = 0.25*(1.0+xi(:,2)).*(2.0*xi(:,1)+xi(:,2));         % dN3(xi,eta)/dxi
dN(:,4) = 0.25*(1.0+xi(:,2)).*(2.0*xi(:,1)-xi(:,2));         % dN4(xi,eta)/dxi
dN(:,5) =-xi(:,1).*(1.0-xi(:,2));                            % dN5(xi,eta)/dxi
dN(:,6) = 0.50*(1.0-xi(:,2).*xi(:,2));                       % dN6(xi,eta)/dxi
dN(:,7) =-xi(:,1).*(1.0+xi(:,2));                            % dN7(xi,eta)/dxi
dN(:,8) =-0.50*(1.0-xi(:,2).*xi(:,2));                       % dN8(xi,eta)/dxi
dN(:,9) = 0.25*(1.0-xi(:,1)).*( xi(:,1)+2.0*xi(:,2));        % dN1(xi,eta)/deta
dN(:,10)= 0.25*(1.0+xi(:,1)).*(-xi(:,1)+2.0*xi(:,2));        % dN2(xi,eta)/deta
dN(:,11)= 0.25*(1.0+xi(:,1)).*( xi(:,1)+2.0*xi(:,2));        % dN3(xi,eta)/deta
dN(:,12)= 0.25*(1.0-xi(:,1)).*(-xi(:,1)+2.0*xi(:,2));        % dN4(xi,eta)/deta
dN(:,13)=-0.50*(1.0-xi(:,1).*xi(:,1));                       % dN5(xi,eta)/deta
dN(:,14)=-xi(:,2).*(1.0+xi(:,1));                            % dN6(xi,eta)/deta
dN(:,15)= 0.50*(1.0-xi(:,1).*xi(:,1));                       % dN7(xi,eta)/deta
dN(:,16)=-xi(:,2).*(1.0-xi(:,1));                            % dN8(xi,eta)/deta


% --- Jacobian ---
% Jac =  [dx/dxi  dx/deta
%         dz/dxi  dz/deta]
Jac=zeros(nXi,2,2);
for iXi=1:nXi
    for iNod=1:8
        Jac(iXi,1,1)=Jac(iXi,1,1)+dN(iXi,iNod)*Node(iNod,1);
        Jac(iXi,1,2)=Jac(iXi,1,2)+dN(iXi,iNod)*Node(iNod,2);
        Jac(iXi,2,1)=Jac(iXi,2,1)+dN(iXi,8+iNod)*Node(iNod,1);
        Jac(iXi,2,2)=Jac(iXi,2,2)+dN(iXi,8+iNod)*Node(iNod,2);
    end
end

% INITIALIZATION
Me=zeros(16,16);
Ke=zeros(16,16);

% --- Integration: loop over integration points ---
for iXi=1:nXi
    
    % INTEGRATION POINT COORDINATES
    xiX=0;
    xiY=0;
    for iNod=1:8
        xiX=xiX+N(iXi,iNod)*Node(iNod,1);
        xiY=xiY+N(iXi,iNod)*Node(iNod,2);
    end
    
    % JACOBIAN
    JacUtil=[Jac(iXi,1,1),Jac(iXi,1,2);Jac(iXi,2,1),Jac(iXi,2,2)];
    detJac=det(JacUtil);
    DN1= JacUtil\[dN(iXi,1);dN(iXi,9)];  % [dN1/dx; dN1/dz]
    DN2= JacUtil\[dN(iXi,2);dN(iXi,10)]; % [dN2/dx; dN2/dz]
    DN3= JacUtil\[dN(iXi,3);dN(iXi,11)]; % [dN3/dx; dN3/dz]
    DN4= JacUtil\[dN(iXi,4);dN(iXi,12)]; % [dN4/dx; dN4/dz]
    DN5= JacUtil\[dN(iXi,5);dN(iXi,13)]; % [dN5/dx; dN5/dz]
    DN6= JacUtil\[dN(iXi,6);dN(iXi,14)]; % [dN6/dx; dN6/dz]
    DN7= JacUtil\[dN(iXi,7);dN(iXi,15)]; % [dN7/dx; dN7/dz]
    DN8= JacUtil\[dN(iXi,8);dN(iXi,16)]; % [dN8/dx; dN8/dz]
    
    if strcmpi(Options.problem,'axisym')
        B=[DN1(1)       0      DN2(1)       0      DN3(1)       0      DN4(1)       0      DN5(1)       0      DN6(1)       0      DN7(1)       0      DN8(1)       0;
           N(iXi,1)/xiX 0      N(iXi,2)/xiX 0      N(iXi,3)/xiX 0      N(iXi,4)/xiX 0      N(iXi,5)/xiX 0      N(iXi,6)/xiX 0      N(iXi,7)/xiX 0      N(iXi,8)/xiX 0;
           0            DN1(2) 0            DN2(2) 0            DN3(2) 0            DN4(2) 0            DN5(2) 0            DN6(2) 0            DN7(2) 0            DN8(2);
           DN1(2)       DN1(1) DN2(2)       DN2(1) DN3(2)       DN3(1) DN4(2)       DN4(1) DN5(2)       DN5(1) DN6(2)       DN6(1) DN7(2)       DN7(1) DN8(2)       DN8(1)];
        Nxi=[N(iXi,1) 0        N(iXi,2) 0        N(iXi,3) 0        N(iXi,4) 0        N(iXi,5) 0        N(iXi,6) 0        N(iXi,7) 0        N(iXi,8) 0;
            0        N(iXi,1) 0        N(iXi,2) 0        N(iXi,3) 0        N(iXi,4) 0        N(iXi,5) 0        N(iXi,6) 0        N(iXi,7) 0        N(iXi,8)];
        
        Me=Me+2*pi*xiX*H(iXi)*detJac*(Nxi.'*rho*Nxi);
        Ke=Ke+2*pi*xiX*H(iXi)*detJac*(B.'*C*B);
    else
        B=[DN1(1)  0       DN2(1)   0       DN3(1)  0      DN4(1)  0      DN5(1)  0      DN6(1)  0      DN7(1)  0      DN8(1)  0;
            0       DN1(2)  0        DN2(2)  0       DN3(2) 0       DN4(2) 0       DN5(2) 0       DN6(2) 0       DN7(2) 0       DN8(2);
            DN1(2)  DN1(1)  DN2(2)   DN2(1)  DN3(2)  DN3(1) DN4(2)  DN4(1) DN5(2)  DN5(1) DN6(2)  DN6(1) DN7(2)  DN7(1) DN8(2)  DN8(1)];
        Nxi=[N(iXi,1) 0        N(iXi,2) 0        N(iXi,3) 0        N(iXi,4) 0        N(iXi,5) 0        N(iXi,6) 0        N(iXi,7) 0        N(iXi,8) 0;
            0        N(iXi,1) 0        N(iXi,2) 0        N(iXi,3) 0        N(iXi,4) 0        N(iXi,5) 0        N(iXi,6) 0        N(iXi,7) 0        N(iXi,8)];
        Me=Me+H(iXi)*detJac*(Nxi.'*rho*Nxi);
        Ke=Ke+H(iXi)*detJac*(B.'*C*B);
    end
end
end
