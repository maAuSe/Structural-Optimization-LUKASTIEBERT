function [F,dFdx] = loads_shell2(DLoad,Node,dDLoaddx,dNodedx)

%LOADS_SHELL2   Equivalent nodal forces for a SHELL2 element in the GCS.
%
%   F = loads_shell2(DLoad,Node)
%   computes the equivalent nodal forces of a distributed load
%   (in the global coordinate system).
%
%   DLoad      Distributed loads     [n1globalX; n1globalY; n1globalZ; ...]
%                                                                   (6 * 1)
%   Node       Node definitions           [x y z] (3 * 3)
%   F          Load vector  (6 * 1)

% Mattias Schevenels
% April 2020

if nargin<3, dDLoaddx = []; end
if nargin<4, dNodedx = []; end
if (nnz(dDLoaddx)+nnz(dNodedx)) > 0
    error('Sensitivities have not been implemented yet.')
end
nVar = 0;
if nargout>1 && (~isempty(dDLoaddx) || ~isempty(dNodedx))
    nVar = max(size(dDLoaddx,3),size(dNodedx,3));
end 

L=norm(Node(2,:)-Node(1,:));

t = trans_shell2(Node);

if size(DLoad,1)>=12 % distributed forces and moments
  T=blkdiag(t,t,t,t);
  DLoadLCS=T*DLoad(1:12,:);
else % only distributed forces
  T=blkdiag(t,t);
  DLoadLCS=T*DLoad(1:6,:);
end

FLCS=loadslcs_shell2(DLoadLCS,L);

T=blkdiag(t,t);
F=T.'*FLCS;

if nargout>1
    dFdx = zeros(size(F,1),size(F,2),nVar); 
end
