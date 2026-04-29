function [Ke,Me,dKedx] = ke_shell8(Node,Section,Material,Options,dNodedx,dSectiondx)

%KE_SHELL8   Shell element stiffness and mass matrix in global coordinate system.
%
%   [Ke,Me] = ke_shell8(Node,Section,Material,Options)
%    Ke     = ke_shell8(Node,Section,Material,Options)
%   returns the element stiffness and mass matrix in the global coordinate system
%   for an eight node shell element.
%
%   Node       Node definitions           [x y z] (8 * 3)
%              Nodes should have the following order:
%              4----7----3
%              |         |
%              8         6
%              |         |
%              1----5----2
%
%   Section    Section definition         [h] or [h1 h2 h3 h4]
%              (uniform thickness or defined in corner nodes(1,2,3,4))
%   Material   Material definition        [E nu rho] or 
%                                         [Exx Eyy nuxy muxy muyz muzx theta rho]
%   Options    Element options struct. Fields:
%              -LCSType: determine the reference local element
%                        coordinate system. Values:
%                        'element' (default) or 'global'
%              -MatType: 'isotropic' (default) or 'orthotropic'
%              -Offset: nodal offset from shell midplane. Values:
%                      'top', 'mid' (default), 'bot' or numerical value
%              -RotaryInertia: 0 (default) or 1
%   Ke         Element stiffness matrix (48 * 48)
%   Me         Element mass matrix (48 * 48)
%
%   This element is based on chapter 15 of
%   The Finite Element Method: for Solid and Structural Mechanics,
%   Zienkiewicz (2005).
%
%   See also KE_BEAM, ASMKM, KE_TRUSS.

% Miche Jansen
% 2009

if nargin<5, dNodedx = []; end
if nargin<6, dSectiondx = []; end
if (nnz(dNodedx)+nnz(dSectiondx)) > 0
    error('Sensitivities have not been implemented yet.')
end
nVar = 0;
if nargout>2 && (~isempty(dNodedx) || ~isempty(dSectiondx))
    nVar = max(size(dNodedx,3),size(dSectiondx,3));
end

if nargout>1, Me = []; end
if nargout>2, dKedx = cell(nVar,1); end


% Check nodes
Node=Node(1:8,1:3);
if ~ all(isfinite(Node))
    error('Not all the nodes exist.')
end


% Section
h=Section(:);
if length(h)<4
    h = h(1);
elseif any(isnan(h(2:4)))
    h = h(1);
else
    h = [h(1:4);(h(1)+h(2))/2;(h(2)+h(3))/2;(h(3)+h(4))/2;(h(1)+h(4))/2];
end
% offset
if isfield(Options,'Offset')
    if isnumeric(Options.Offset)
        d = Options.Offset;
    else
        switch lower(Options.Offset)
            case 'top'; d = h/2;
            case 'mid'; d = 0;
            case 'bot'; d = -h/2;
        end
    end
else
    d = 0;
end
if isfield(Options,'RotaryInertia')
    RotaryInertia = Options.RotaryInertia;
else
    RotaryInertia = 0;
end

% integratiepunten
[x,H] = gaussq(2);
[xzeta,Hzeta] = gaussq(2);


%bepaling k-factor (formule op p14-222 van ansys theorie)
% A=0;
% for iGauss=1:size(x,1)
%     xi = x(iGauss,1);
%     eta = x(iGauss,2);
%     [Ni,dN_dxi,dN_deta] = sh_qs8(xi,eta);
% f = norm([Node(:,2).'*dN_dxi*Node(:,3).'*dN_deta-Node(:,2).'*dN_deta*Node(:,3).'*dN_dxi, ...
%      Node(:,3).'*dN_dxi*Node(:,1).'*dN_deta-Node(:,3).'*dN_deta*Node(:,1).'*dN_dxi, ...
%      Node(:,1).'*dN_dxi*Node(:,2).'*dN_deta-Node(:,1).'*dN_deta*Node(:,2).'*dN_dxi]);
%  A = A + f;
% end
A = norm(cross(Node(2,:)-Node(1,:),Node(4,:)-Node(1,:)));
k=max([1.2 1+0.2*A/(25*mean(h)^2)]);
% k = 1.2;

% Material
E=Material(1,1);
if isfield(Options,'MatType')
    MatType = Options.MatType;
else
    MatType = 'isotropic';
end
if nargout==2
    if strcmpi(MatType,'isotropic')
        rho=Material(1,3);
    else rho=Material(1,8);
    end
end
D = cmat_shell8(MatType,Material,k);

xi_eta = [ -1 -1;1  -1;1 1; -1 1;0 -1; 1 0; 0 1; -1 0];
v3i = zeros(8,3);
for iNode = 1:8
    [Ni,dN_dxi,dN_deta] = sh_qs8(xi_eta(iNode,1),xi_eta(iNode,2));
    Jm = [dN_dxi.'*Node(:,1) dN_dxi.'*Node(:,2) dN_dxi.'*Node(:,3);
        dN_deta.'*Node(:,1) dN_deta.'*Node(:,2) dN_deta.'*Node(:,3)];
    v3i(iNode,:) = cross(Jm(1,:),Jm(2,:));
    v3i(iNode,:) = v3i(iNode,:)/norm(v3i(iNode,:));
end

% referentie lokaal assenkruis
v0 = trans_shell8(Node,Options);

% lokaal assenkruis per knoop --> projectie-methode COOK p579
v1i = zeros(8,3);
v2i = zeros(8,3);
for iNode=1:8
    v3i(iNode,:) = v3i(iNode,:)/norm(v3i(iNode,:));
    v2i(iNode,:) = cross(v3i(iNode,:),v0(1,:));
    if norm(v2i(iNode,:))<=0.02
        v1i(iNode,:)= cross(v0(2,:),v3i(iNode,:));
        v2i(iNode,:)= cross(v3i(iNode,:),v1i(iNode,:));
    else
        v1i(iNode,:)= cross(v2i(iNode,:),v3i(iNode,:));
    end
    v1i(iNode,:) = v1i(iNode,:)/norm(v1i(iNode,:));
    v2i(iNode,:) = v2i(iNode,:)/norm(v2i(iNode,:));
end

% print nodal coordinate systems
% fprintf('\n%g %g %g %g %g %g %g %g \n %g %g %g %g %g %g %g %g \n %g %g %g %g %g %g %g %g \n',v1i)

Ke = zeros(48,48);
if nargout==2
    zetan = 2*d./h;
    Me = zeros(48,48);
end

Aelem = 0;

for iGauss=1:size(x,1)
    xi = x(iGauss,1);
    eta = x(iGauss,2);
    
    % 2d vormfuncties
    [Ni,dN_dxi,dN_deta] = sh_qs8(xi,eta);
    
    %bepaling lokale z-richting bij het gausspunt
    %     Jm = [dN_dxi.'*Node(:,1) dN_dxi.'*Node(:,2) dN_dxi.'*Node(:,3);
    %          dN_deta.'*Node(:,1) dN_deta.'*Node(:,2) dN_deta.'*Node(:,3)];
    %     v3g=cross(Jm(1,:),Jm(2,:));
    %     [mini,Ig]=min(v3g);
    %     v1g = cross(Ei(Ig,:),v3g);
    %     v2g = cross(v3g,v1g);
    v3g = Ni.'*v3i;
    v2g = Ni.'*v2i;
    v1g = Ni.'*v1i;
    
    %transformatie matrix lokaal <-> globaal
    A = [v1g/norm(v1g); v2g/norm(v2g); v3g/norm(v3g)];
    
    % transformatie matrix uit Cook p274 en p194 Zienkiewicz deel1
    theta = vtrans_solid(A,'strain');
    theta = theta([1,2,4,5,6],:);
    %materiaalmatrix in globale assenkruis (p198 Zienkiewicz deel1)
    Dg = theta.'*D*theta;
    
    for zeta = [-xzeta(1),xzeta(1)]
        
        [Bg,J] = b_shell8(Ni,dN_dxi,dN_deta,zeta,Node,h,v1i,v2i,v3i,d);
        %B = theta*Bg;
        
        Ke = Ke + Bg.'*Dg*Bg*H(iGauss)*det(J);
        
        if nargout==2
            N = zeros(3,48);
            N(1,1:6:43) = Ni;
            N(2,2:6:44) = Ni;
            N(3,3:6:45) = Ni;
            if RotaryInertia == 1  
                N(:,4:6:46) = (Ni(:,[1 1 1]).*(diag((zeta-zetan).*(h/2))*(-v2i))).';
                N(:,5:6:47) = (Ni(:,[1 1 1]).*(diag((zeta-zetan).*(h/2))*(v1i))).';
            end
            Me = Me +rho*det(J)*(N.'*N)*H(iGauss);
        end
    end
    
    Aelem2 = norm([Node(:,2).'*dN_dxi*Node(:,3).'*dN_deta-Node(:,2).'*dN_deta*Node(:,3).'*dN_dxi, ...
        Node(:,3).'*dN_dxi*Node(:,1).'*dN_deta-Node(:,3).'*dN_deta*Node(:,1).'*dN_dxi, ...
        Node(:,1).'*dN_dxi*Node(:,2).'*dN_deta-Node(:,1).'*dN_deta*Node(:,2).'*dN_dxi]);
    Aelem = Aelem + Aelem2;
end

% oplossing drilling COOK p575
Ke(6:6:48,6:6:48)=(1e-5)*E*Aelem*mean(h)*((8/7)*eye(8)-(1/7)*ones(8,8));

% % eenvoudige oplossing voor drilling uit BATHE p209
% KeD = diag(Ke);
% minKeD = min(KeD(KeD ~= 0));
% for e=6:6:48
% Ke(e,e) = minKeD/1000;
% end


%rotatiecomponenten transformeren naar globaal assenstelsel
T = eye(48);
for iNode=1:8
T((4:6)+6*(iNode-1),(4:6)+6*(iNode-1)) = [ v2i(iNode,:); v1i(iNode,:); v3i(iNode,:)];
end

Ke = T.'*Ke*T;
if nargout==2
    Me = T.'*Me*T';
end

end

function z = cross(x,y)
z = x([2 1 1]).*y([3 3 2])-y([2 1 1]).*x([3 3 2]);
z(2) = -z(2);
end