function [Ke,Me,dKedx] = ke_shell2(Node,Section,Material,Options,dNodedx,dSectiondx)

%KE_SHELL2   SHELL2 element stiffness matrix in global coordinate system.
%
%   Ke = ke_shell2(Node,Section,Material) returns the element stiffness
%   matrix in the global coordinate system for a two-node axisymmetric shell
%   element. The global y-axis is assumed to be the axis of symmetry.
%
%   Node       Node definitions           [x y z] (3 * 3)
%   Section    Section definition         [h]
%   Material   Material definition        [E nu]
%   Ke         Element stiffness matrix (6 * 6)

% Mattias Schevenels
% April 2020

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

% CHECK NODES
if ~all(isfinite(Node(1:2,1:3)))
  error('Not all the nodes exist.');
end
if any(Node(:,3)~=0) || any(Node(:,1)<0)
  error('Axisymmetric SHELL2 elements must be defined in the x,y-plane with x>=0.');
end
if all(Node(:,1)==0)
  error('Axisymmetric SHELL2 elements must no be defined along the y-axis.');
end

% ELEMENT LENGTH, LOCATION, SLOPE
L = norm(Node(2,:)-Node(1,:));
r1 = Node(1,1);
phi = atan2(Node(2,2)-Node(1,2),Node(2,1)-Node(1,1));

% MATERIAL
E = Material(1);
nu = Material(2);

% SECTION
h = Section(1);

% TRANSFORMATION MATRIX
t = trans_shell2(Node);
T = blkdiag(t,t);

% STIFFNESS MATRIX
    KeLCS = kelcs_shell2(r1,phi,L,h,E,nu);
    Ke = T.'*KeLCS*T;
end
