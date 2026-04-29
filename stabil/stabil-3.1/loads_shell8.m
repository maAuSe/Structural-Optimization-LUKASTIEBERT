function [F,dFdx] = loads_shell8(DLoad,Node,dDLoaddx,dNodedx)

%LOADS_SHELL8   Equivalent nodal forces for a shell8 element in the GCS.
%
%   F = loads_shell8(DLoad,Node)
%   computes the equivalent nodal forces of a distributed load 
%   (in the global coordinate system).
%
%   DLoad      Distributed loads      [n1globalX; n1globalY; n1globalZ; ...]
%              in corner Nodes         (12 * 1)
%   Node       Node definitions       [x y z] (8 * 3)
%   F          Load vector  (48 * 1)
%
%   See also LOADSLCS_BEAM, ELEMLOADS, LOADS_TRUSS.

% Miche Jansen
% 2009

if nargin<3, dDLoaddx = []; end
if nargin<4, dNodedx = []; end
if (nnz(dDLoaddx)+nnz(dNodedx)) > 0
    error('Sensitivities have not been implemented yet.')
end
nVar = 0;
if nargout>1 && (~isempty(dDLoaddx) || ~isempty(dNodedx))
    nVar = max(size(dDLoaddx,3),size(dNodedx,3));
end 


Node = Node(1:8,1:3);
DLoad = DLoad(1:12,:);

% integratiepunten
[x,H] = gaussq(3); %moet drie zijn ingeval de belasting varieert over het element?

F = zeros(48,size(DLoad,2));
DLoad = [DLoad;(DLoad((1:3),:)+DLoad((4:6),:))/2;(DLoad((4:6),:)+DLoad((7:9),:))/2;...
    (DLoad((7:9),:)+DLoad((10:12),:))/2;(DLoad((10:12),:)+DLoad((1:3),:))/2];

%p148 Zienkiewicz deel1 (2005)
for iGauss=1:size(x,1)
    xi = x(iGauss,1);
    eta = x(iGauss,2);
    [Ni,dN_dxi,dN_deta] = sh_qs8(xi,eta);
    N = zeros(3,48);
    N(1,1:6:43) = Ni;
    N(2,2:6:44) = Ni;
    N(3,3:6:45) = Ni;
    Dloadg = [Ni.'*DLoad(1:3:22,:); Ni.'*DLoad(2:3:23,:); Ni.'*DLoad(3:3:24,:)];
    det = norm(cross(Node.'*dN_dxi,Node.'*dN_deta));

    F = F + H(iGauss)*N.'*Dloadg*det;
end

if nargout>1
    dFdx = zeros(size(F,1),size(F,2),nVar);  % avoid unassignment error
end
