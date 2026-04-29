function SeLCS = selcs_plane6(Node,Section,Material,UeLCS,Options)

%SELCS_PLANE6   Compute the element stresses for a plane6 element.
%
%   [SeLCS] = selcs_plane6(Node,Section,Material,UeGCS,Options)
%   computes the element stresses in the
%   local coordinate system for the plane6 element.
%
%   Node       Node definitions           [x y z] (6 * 3)
%              Nodes should have the following order:
%              3
%              | \
%              6  5
%              |    \
%              1--4--2
%   Section    Section definition  [h]  (only used in plane stress)
%   Material   Material definition [E nu rho]
%   UeLCS      Displacements (6 * nSteps)
%   Options    Element options            {Option1 Option2 ...}
%   SeLCS      Element stresses in LCS in corner nodes IJKL
%              (9 * nTimeSteps) [sxx syy sxy]
%
%   See also ELEMSTRESS, SE_PLANE6.

% Jacinto Ulloa
% 2017

% Check nodes
if ~all(isfinite(Node(1:6,1:3))), error('Not all the nodes exist.'); end

% Material
E = Material(1,1);
nu = Material(1,2);

% Options default values
if nargin < 3, Options = []; end
if ~isfield(Options,'problem'),      Options.problem = 'stress'; end
if ~isfield(Options,'bendingmodes'), Options.bendingmodes = true; end

% Lame parameters
lambda = nu*E/(1+nu)/(1-2*nu);
mu = E/2/(1+nu);

% Constitutive matrix
if any(strcmpi(Options,'axisym'))           % Axisymmetric
  D=[lambda+2*mu lambda      lambda       0
     lambda      lambda+2*mu lambda       0
     lambda      lambda      lambda+2*mu  0
     0           0           0            mu];
elseif any(strcmpi(Options,'2dstrain'))  % Plane strain
  D=[lambda+2*mu lambda      0   
     lambda      lambda+2*mu 0   
     0           0           mu];
else                                        % Plane stress
  D=[ 4*mu-4*mu^2/(lambda+2*mu)   2*lambda*mu/(lambda+2*mu)     0
      2*lambda*mu/(lambda+2*mu)   4*mu-4*mu^2/(lambda+2*mu)     0
                             0                           0     mu];
end

% Compute the stress directly at the corner nodes
% Corner nodes in natural coordinates
% xi(1,1) = 0; eta(1,1) = 0;
% xi(1,2) = 1; eta(1,2) = 0;
% xi(1,3) = 0; eta(1,3) = 1;
% SeLCS = zeros(1,9); 
% for iNod = 1:3
%         [~,dNi_dxi,dNi_deta]=sh_t6(xi(iNod),eta(iNod));
%         J = [dNi_dxi dNi_deta]'*Node(1:6,1:2);
%         dN = J\[dNi_dxi dNi_deta]';
%         B = zeros(3,12);
%         B(1,1:2:11) = dN(1,:);
%         B(2,2:2:12) = dN(2,:);
%         B(3,1:2:11) = dN(2,:);
%         B(3,2:2:12) = dN(1,:);
%         SeLCS(1,(3*iNod-2):3*iNod) = (D*B*UeLCS);
% end
% SeLCS = (SeLCS)';

% Compute the stress at the Gauss points and extrapolate to the corners
% Gauss points in natural coordinates
xig(1) = 0.16666666666667; etag(1) = 0.16666666666667;
xig(2) = 0.66666666666667; etag(2) = 0.16666666666667;
xig(3) = 0.16666666666667; etag(3) = 0.66666666666667;
SeGP = zeros(1,9); 
for iNod = 1:3
        [~,dNi_dxi,dNi_deta]=sh_t6(xig(iNod),etag(iNod));
        J = [dNi_dxi dNi_deta]'*Node(1:6,1:2);
        dN = J\[dNi_dxi dNi_deta]';
        B = zeros(3,12);
        B(1,1:2:11) = dN(1,:);
        B(2,2:2:12) = dN(2,:);
        B(3,1:2:11) = dN(2,:);
        B(3,2:2:12) = dN(1,:);
        SeGP(1,(3*iNod-2):3*iNod) = (D*B*UeLCS);
end
SeGP = (SeGP)';
% Corner nodes in the xi', eta' coordinate system
xic(1) = -0.33333333333333; etac(1) = -0.33333333333333;
xic(2) = 1.66666666666667; etac(2) = -0.33333333333333;
xic(3) = -0.33333333333333; etac(3) = 1.66666666666667;
% Matrix of shape functions evaluated at the corner points
shape = [1-xic(1,1)-etac(1,1) xic(1) etac(1);
         1-xic(1,2)-etac(2) xic(2) etac(2);
         1-xic(1,3)-etac(3) xic(3) etac(3)];
SeLCS = zeros(9,1);
for iDir = 1:3
    SeLCS([iDir iDir+3 iDir+6],1) = shape*SeGP([iDir iDir+3 iDir+6],1);
end
end

