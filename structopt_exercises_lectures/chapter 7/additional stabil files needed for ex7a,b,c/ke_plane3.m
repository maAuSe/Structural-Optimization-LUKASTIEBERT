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
XY = [1 X(1) Y(1);
      1 X(2) Y(2);
      1 X(3) Y(3)];

Delta=0.5*det(XY);

% Shape function derivatives
BC = [b1  0  b2   0  b3   0
      0   c1  0  c2   0  c3
      c1  b1 c2  b2  c3  b3];
Be=1/(2*Delta)*BC;

% Stiffness matrix
Ke=Delta*Be.'*C*Be;

% Stiffness matrix sensitivities
if nargout>=3 && (any(dNodedx(:)~=0) || any(dSectiondx(:)~=0))
  nVar = 0;
  if nargout>2 && (~isempty(dNodedx) || ~isempty(dSectiondx))
    nVar = max(size(dNodedx,3),size(dSectiondx,3));
  end
  if nVar==0 || isempty(dNodedx), dNodedx = zeros([size(Node),nVar]); end
  if nVar==0 || isempty(dSectiondx), dSectiondx = zeros([size(Section),nVar]); end

  if strcmpi(Options.problem,'2dstress')
    dCdx = C/Section.*dSectiondx;
  else
    dCdx = zeros([size(C),nVar]);
  end

  dXdx = dNodedx(:,1,:);
  dYdx = dNodedx(:,2,:);
  db1dx = dYdx(2,1,:)-dYdx(3,1,:); db2dx = dYdx(3,1,:)-dYdx(1,1,:); db3dx = dYdx(1,1,:)-dYdx(2,1,:);
  dc1dx = dXdx(3,1,:)-dXdx(2,1,:); dc2dx = dXdx(1,1,:)-dXdx(3,1,:); dc3dx = dXdx(2,1,:)-dXdx(1,1,:);

  nul = zeros(1,1,nVar);
  dXYdx = [nul dXdx(1,1,:) dYdx(1,1,:);
           nul dXdx(2,1,:) dYdx(2,1,:);
           nul dXdx(3,1,:) dYdx(3,1,:)];
  dDeltadx = zeros(1,1,nVar);
  for iVar = 1:nVar
    dDeltadx(iVar) = 0.5*trace(adjoint(XY)*dXYdx(:,:,iVar));
  end
  dBCdx = [db1dx  nul   db2dx  nul    db3dx  nul
           nul    dc1dx nul    dc2dx  nul    dc3dx
           dc1dx  db1dx dc2dx  db2dx  dc3dx  db3dx];
  dBedx = -0.5*dDeltadx/Delta^2.*BC+0.5/Delta*dBCdx;
  dKedx = cell(nVar,1);
  for iVar = 1:nVar
    dKedx{iVar} = dDeltadx(:,:,iVar)*Be.'*C*Be + Delta*dBedx(:,:,iVar).'*C*Be + Delta*Be.'*dCdx(:,:,iVar)*Be + Delta*Be.'*C*dBedx(:,:,iVar);
  end
else
  dKedx = repmat({zeros(6)},nVar,1);
end

% Mass matrix
if (nargout==2)
    Me=1/3*rho*Delta*eye(6);
end
end
