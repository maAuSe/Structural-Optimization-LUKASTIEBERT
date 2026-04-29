function [Ke,Me,dKedx] = ke_mass(Node,Section,Material,Options,dNodedx,dSectiondx)
%KE_MASS   mass element system matrices in the global coordinate system.
%
%   [Ke,Me] = ke_mass(Node,Section,Material,Options) returns the element
%   stiffness and mass matrix in the global coordinate system for a 
%   concentrated mass element
%
%   Node       Node definitions           [x1 y1 z1] (1 * 3)
%   Section    Section definition         [m]
%   Material   Material definition        []
%   Options    Element options            {Option1 Option2 ...}
%   Ke         Element stiffness matrix   (6 * 6)
%   Me         Element mass matrix        (6 * 6)

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

% Check nodes
if ~ all(isfinite(Node(1,1:3)))
    error('Not all the nodes exist.')
end

% Section
m=Section(1);

Ke=zeros(6);
Me=diag([m m m 0 0 0]);
end
