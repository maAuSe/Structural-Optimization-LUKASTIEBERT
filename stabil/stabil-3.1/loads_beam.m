function [F,dFdx] = loads_beam(DLoad,Node,dDLoaddx,dNodedx)

%LOADS_BEAM   Equivalent nodal forces for a beam element in the GCS.
%
%   F = LOADS_BEAM(DLoad,Node)
%   computes the equivalent nodal forces of a distributed load 
%   (in the global coordinate system).
%
%   [F,dFdx] = LOADS_BEAM(DLoad,Node,dDLoaddx,dNodedx)
%   additionally computes the derivatives of the equivalent nodal forces
%   with respect to the design variables x.
%
%   DLoad      Distributed loads    [n1globalX; n1globalY; n1globalZ; ...]  
%                   (6/8 * nLC * nDLoads)
%   Node       Node definitions     [x y z] (3 * 3)
%   dDLoaddx   Distributed loads derivatives    (6/8 * nLC * nVar * nDLoads)
%   dNodedx    Node definitions derivatives     (SIZE(Node) * nVar)
%   F          Load vector              (12 * nLC)
%   dFdx       Load vector derivatives  (12 * nLC * nVar)
%
%   See also LOADSLCS_BEAM, ELEMLOADS, LOADS_TRUSS.

% David Dooms, Wouter Dillen
% March 2008, April 2017

if nargin<3, dDLoaddx = []; end
if nargin<4, dNodedx = []; end

nVar = 0;
if nargout>1 && (~isempty(dDLoaddx) || ~isempty(dNodedx))
    nVar = max(size(dDLoaddx,3),size(dNodedx,3));
end 

if nVar==0 || isempty(dDLoaddx), dDLoaddx = zeros(size(DLoad,1),size(DLoad,2),nVar,size(DLoad,3)); end
if nVar==0 || isempty(dNodedx), dNodedx = zeros([size(Node),nVar]); end


% Element length
L=norm(Node(2,:)-Node(1,:));
dLdx=(permute(dNodedx(2,:,:)-dNodedx(1,:,:),[3 2 1])*(Node(2,:)-Node(1,:)).')/L;

% transformation matrix
[t,dtdx] = trans_beam(Node,dNodedx);
T=blkdiag(t,t);
dTdx = [dtdx  zeros(size(dtdx));
            zeros(size(dtdx)) dtdx];

% DLoadLCS
[DLoadLCS,dDLoadLCSdx] = dloadgcs2lcs(T,DLoad,dTdx,dDLoaddx);
    
% FLCS
FLCS=zeros(12,size(DLoad,2));
dFLCSdx=zeros(size(FLCS,1),size(FLCS,2),nVar);
nDLoadsOnElement = size(DLoad,3);
for i=1:nDLoadsOnElement
    [tmpFLCS,tmpdFLCSdx] = loadslcs_beam(DLoadLCS(:,:,i),L,dDLoadLCSdx(:,:,:,i),dLdx);
    FLCS = FLCS + tmpFLCS;
    dFLCSdx = dFLCSdx + tmpdFLCSdx;
end

% F
T=blkdiag(t,t,t,t);
dTdx = [dTdx  zeros(size(dTdx));
        zeros(size(dTdx)) dTdx];
F=T.'*FLCS; 
dFdx=zeros(size(F,1),size(F,2),nVar);
for n=1:nVar
    if any(reshape(dtdx(:,:,n),[],1)~=0)
        dFdx(:,:,n) = dTdx(:,:,n).'*FLCS + T.'*dFLCSdx(:,:,n);
    else
        dFdx(:,:,n) = T.'*dFLCSdx(:,:,n);
    end
end




