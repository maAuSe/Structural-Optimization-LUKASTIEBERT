function [F,dFdx] = loads_shell4(DLoad,Node,dDLoaddx,dNodedx)

%LOADS_SHELL4   Equivalent nodal forces for a shell4 element in the GCS.
%
%   F = loads_shell4(DLoad,Node)
%   computes the equivalent nodal forces of a distributed load 
%   (in the global coordinate system).
%
%   DLoad      Distributed loads      [n1globalX; n1globalY; n1globalZ; ...]
%              in corner Nodes         (12 * 1)
%   Node       Node definitions       [x y z] (4 * 3)
%   F          Load vector  (24 * 1)
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


Node = Node(1:4,1:3);
DLoad = DLoad(1:12,:);

[t,Node_lc]=trans_shell4(Node);

T=blkdiag(t,t,t,t);

DLoadLCS=T*DLoad;

FLCS=loadslcs_shell4(DLoadLCS,Node_lc);

T=blkdiag(t,t,t,t,t,t,t,t);

F=T.'*FLCS; 


if nargout>1
    dFdx = zeros(size(F,1),size(F,2),nVar);  % avoid unassignment error
end
end