function [Ke,Me,dKedx] = ke_shell4(Node,Section,Material,Options,dNodedx,dSectiondx)

%KE_SHELL4   shell element stiffness and mass matrix in global coordinate system.
%
%   [Ke,Me] = ke_shell4(Node,Section,Material,Options) 
%    Ke     = ke_shell4(Node,Section,Material,Options) 
%   returns the element stiffness and mass matrix in the global coordinate system
%   for a four node shell element (isotropic material).
%
%   Node       Node definitions           [x y z] (4 * 3)
%              Nodes should have the following order:
%              4---------3
%              |         |
%              |         |
%              |         |
%              1---------2
%
%   Section    Section definition         [h]
%   Material   Material definition        [E nu rho]
%   Options    Element options            {Option1 Option2 ...}
%   Ke         Element stiffness matrix (24 * 24)
%   Me         Element mass matrix (24 * 24)
%
%   This shell element consists of a bilinear membrane element and 
%   four overlaid DKT triangles for the bending stiffness.
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
Node=Node(1:4,1:3);
if ~ all(isfinite(Node))
    error('Not all the nodes exist.')
end

% Material
E=Material(1,1);
nu=Material(1,2);

% Section
h=Section(1,1);

% Transformation matrix
[t,Node_lc,W]=trans_shell4(Node);
T=blkdiag(t,t,t,t,t,t,t,t);


if nargout==2     % stiffness and mass
    if nargin<4
        Options={};
    end
    rho=Material(1,3);
    [KeLCS,MeLCS]=kelcs_shell4(Node_lc,h,E,nu,rho,Options);
    Me=T.'*W*MeLCS*W.'*T;
    Ke=T.'*W*KeLCS*W.'*T;
else                            % only stiffness
    KeLCS=kelcs_shell4(Node_lc,h,E,nu);
    Ke=T.'*W*KeLCS*W.'*T;
end



end
