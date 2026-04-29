function [Forces,dForcesdx]=forceslcs_truss(KeLCS,UeLCS,dTm,A,E,alpha,dKeLCSdx,dUeLCSdx)

%FORCESLCS_TRUSS   Compute the element forces for a truss element in the LCS.
%
%   Forces = FORCESLCS_TRUSS(KeLCS,UeLCS,dTm,A,E,alpha)
%   Forces = FORCESLCS_TRUSS(KeLCS,UeLCS)
%       computes the element forces for the truss element in the local coordinate 
%       system (algebraic convention).
%
%   [Forces,dForcesdx] = FORCESLCS_TRUSS(KeLCS,UeLCS,[],[],[],[],dKeLCSdx,dUeLCSdx)
%       additionally computes the derivatives of the element forces with
%       respect to the design variables x.
%
%   KeLCS      Element stiffness matrix (6 * 6)
%   UeLCS      Displacements            (6 * nLC)
%   dKeLCSdx   Element stiffness matrix derivatives    (CELL(nVar,1))
%   dUeLCSdx   Displacements derivatives               (SIZE(UeLCS) * nVar)
%   Forces     Element forces                          [N; 0; 0](6 * nLC)
%   dForcesdx  Element forces derivatives              (6 * nLC * nVar)
%
%   See also FORCES_TRUSS, FORCESLCS_BEAM, ELEMFORCES

% David Dooms, Wouter Dillen
% October 2008, April 2017

if nargin<3, dTm = []; end
if nargin<4, A = []; end
if nargin<5, E = []; end
if nargin<6, alpha = []; end
if nargin<7, dKeLCSdx = []; end
if nargin<8, dUeLCSdx = []; end

nVar = 0;
if nargout>1 && (~isempty(dKeLCSdx) || ~isempty(dUeLCSdx))
    nVar = max([length(dKeLCSdx),size(dUeLCSdx,3)]);
end

if nVar==0 || isempty(dKeLCSdx), dKeLCSdx = cell(nVar,1); dKeLCSdx(:) = {zeros(size(KeLCS))}; end
if nVar==0 || isempty(dUeLCSdx), dUeLCSdx = zeros(size(UeLCS,1),size(UeLCS,2),nVar); end

dForcesdx=zeros(size(dUeLCSdx)); 


% compute element forces
Forces=KeLCS*UeLCS;
for n=1:nVar
    dForcesdx(:,:,n)=dKeLCSdx{n}*UeLCS+KeLCS*dUeLCSdx(:,:,n);
end

if ~isempty(dTm)
    % subtract equivalent nodal forces (temperature loads) from element forces
    Flcs=zeros(6,size(dTm,2));
    Flcs(1,:)=-E*A*alpha*dTm;
    Flcs(4,:)= E*A*alpha*dTm;
    Forces=Forces-Flcs;
end

