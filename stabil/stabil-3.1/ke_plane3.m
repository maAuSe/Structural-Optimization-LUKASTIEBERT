function [Ke,Me,dKedx]=ke_plane3(Node,Section,Material,Options,dNodedx,dSectiondx)

%KE_PLANE3   plane element stiffness and mass matrix in global coordinate system.
%
%   [Ke,Me] = ke_plane3(Node,Section,Material,Options) returns the element
%   stiffness and mass matrix in the global coordinate system for a 3-node
%   CST element. Plane3 only operates in the 2D xy-plane so that
%   z-coordinates should be equal to zero.
%
%   Node       Node definitions           [x y z] (3 * 3)
%   Section    Section definitions        [h] (only used in plane stress)
%   Material   Material definition        [E nu rho]
%   Options    Struct containing optional parameters. Fields:
%      .problem Plane stress or plane strain
%               {'2dstress' (default) | '2dstrain'}
%   Ke         Element stiffness matrix (6 * 6)
%   Me         Element mass matrix (6 * 6)

%   Stijn François
%   2016

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

if nargout>1, Me = []; end
if (nargout==2), rho=Material(3); end

% Constitutive matrix
C=cmat_isotropic(Options.problem,Section,Material);

% Triangle shape function
X=Node(:,1);
Y=Node(:,2);
b1=Y(2)-Y(3); b2=Y(3)-Y(1); b3=Y(1)-Y(2);
c1=X(3)-X(2); c2=X(1)-X(3); c3=X(2)-X(1);

% Element area
Delta=0.5*det([1 X(1) Y(1) 
               1 X(2) Y(2) 
               1 X(3) Y(3)]);

% Shape function derivatives
Be=1/(2*Delta)*[b1  0  b2   0  b3   0 
                0   c1  0  c2   0  c3
                c1  b1 c2  b2  c3  b3];

% Stiffness matrix
Ke=Delta*Be.'*C*Be;

% Mass matrix
if (nargout==2)
    Me=1/3*rho*Delta*eye(6);
end
end
