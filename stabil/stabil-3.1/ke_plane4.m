function [Ke,Me,dKedx]=ke_plane4(Node,Section,Material,Options,dNodedx,dSectiondx)

%KE_PLANE4   plane element stiffness and mass matrix in global coordinate system.
%
%   [Ke,Me] = ke_plane4(Node,Section,Material,Options) returns the element
%   stiffness and mass matrix in the global coordinate system for a 4-node
%   plane element. Plane4 only operates in the 2D xy-plane so that
%   z-coordinates should be equal to zero.
%
%   Node       Node definitions           [x y z] (4 * 3)
%              Nodes should have the following order:
%
%                 4---3
%                 |   |
%                 1---2
%
%   Section    Section definition         [h] (only used in plane stress)
%   Material   Material definition        [E nu rho]
%   Options    Struct containing optional parameters. Fields:
%      .problem      Plane stress, plane strain or axisymmetrical
%                    {'2dstress' (default) | '2dstrain' | 'axisym'}
%      .bendingmodes Include (non-conforming) bending modes
%                    {true (default) | false}
%      .integration  Full (2x2) or reduced (1x1) integration
%                    {'full' (default) | 'reduced'}
%   Ke         Element stiffness matrix (8 * 8)
%   Me         Element mass matrix (8 * 8)
%
%   See also KE_BEAM, ASMKM, KE_TRUSS.

% Miche Jansen, Stijn François
% 2016

% Options
if nargin < 4, Options = []; end
if ~isfield(Options,'problem'),      Options.problem = '2dstress'; end
if ~isfield(Options,'bendingmodes'), Options.bendingmodes = ~strcmpi(Options.problem,'axisym'); end
if ~isfield(Options,'integration'),  Options.integration = 'full'; end

if nargin < 5, dNodedx = []; end
if nargin < 6, dSectiondx = []; end
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

if nargout>1, Me = []; end
if (nargout==2), rho=Material(3); end

K11 = zeros(8,8);
if Options.bendingmodes
    K21=zeros(4,8);
    K22=zeros(4,4);
    K12=zeros(8,4);
    J0=[0.25*[-1 1 1 -1]*Node(:,[1,2]);0.25*[-1 -1 1 1]*Node(:,[1,2])];
    detJ0=det(J0);
end

Ke=zeros(8,8);
if (nargout==2), Me=zeros(8,8); end

switch lower(Options.integration)
    case 'reduced'
        xi=[0 0];
        H=[4];
    otherwise % Full integration: the default
        [xi,H]=gaussq(2);
end
nXi=numel(H);

C=cmat_isotropic(Options.problem,Section,Material);

switch lower(Options.problem)
    case {'2dstress', '2dstrain'}
        for iXi=1:nXi
            [Ni,dN_dxi,dN_deta]=sh_qs4(xi(iXi,1),xi(iXi,2));
            J=[dN_dxi.' *Node(:,1)  dN_dxi.' *Node(:,2);
               dN_deta.'*Node(:,1)  dN_deta.'*Node(:,2)];
            detJ = det(J);
            dNi = J\[dN_dxi dN_deta].';
            
            B = zeros(3,8);
            B(1,1:2:7) = dNi(1,:);
            B(2,2:2:8) = dNi(2,:);
            B(3,1:2:7) = dNi(2,:);
            B(3,2:2:8) = dNi(1,:);
            
            K11 = K11 + B.'*C*B*H(1,iXi)*detJ;
            
            % Bending mode contribution
            if Options.bendingmodes
                Nex = [1-xi(iXi,1)^2 -2*xi(iXi,1) 0;
                       1-xi(iXi,2)^2 0 -2*xi(iXi,2)]; % [N5 N6]
                dNex = J0\((detJ0/detJ)*Nex(:,2:3).');
                Bex = zeros(3,4);
                Bex(1,[1 3])=dNex(1,:);
                Bex(2,[2 4])=dNex(2,:);
                Bex(3,[1 3])=dNex(2,:);
                Bex(3,[2 4])=dNex(1,:);
                K21 = K21 + Bex.'*C*B  *H(1,iXi)*detJ;
                K22 = K22 + Bex.'*C*Bex*H(1,iXi)*detJ;
                K12 = K12 + B.'  *C*Bex*H(1,iXi)*detJ;
            end
            
            if (nargout==2)
                N = zeros(2,8);
                N(1,1:2:7) = Ni;
                N(2,2:2:8) = Ni;
                Me = Me + N.'*rho*N*H(1,iXi)*detJ;
            end
        end
        
        if Options.bendingmodes
            Ke = K11-K12*(K22\K21);
        else
            Ke = K11;
        end
        
    case 'axisym'
        for iXi=1:nXi
            [Ni,dN_dxi,dN_deta]=sh_qs4(xi(iXi,1),xi(iXi,2));
            
            J=[dN_dxi.' *Node(:,1) dN_dxi.' *Node(:,2);
               dN_deta.'*Node(:,1) dN_deta.'*Node(:,2)];
            
            dNi = J\[dN_dxi dN_deta].';
            detJ = det(J);
            
            r = Ni.'*Node(:,1);
            
            B = zeros(4,8);
            B(1,1:2:7) = dNi(1,:);   %epsilon_r
            B(2,1:2:7) = Ni./r;      %epsilon_theta
            B(3,2:2:8) = dNi(2,:);   %epsilon_z
            B(4,1:2:7) = dNi(2,:);   %gamma_rz
            B(4,2:2:8) = dNi(1,:);   %gamma_rz
            
            K11 = K11 + 2*pi*r*B.'*C*B*H(1,iXi)*detJ;
            
            if Options.bendingmodes
                warning('No bending modes used in the axisymmetric case.');
            end
            
            if (nargout==2)
                N = zeros(2,8);
                N(1,1:2:7) = Ni;
                N(2,2:2:8) = Ni;
                Me = Me+2*pi*r*N.'*rho*N*H(1,iXi)*detJ;
            end
        end
        Ke = K11;
end
end
