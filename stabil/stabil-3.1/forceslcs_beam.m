function [Forces,dForcesdx] = forceslcs_beam(KeLCS,UeLCS,DLoadLCS,L,TLoadLCS,A,E,alpha,Iyy,Izz,hy,hz,dKeLCSdx,dUeLCSdx,dDLoadLCSdx,dLdx)

%FORCESLCS_BEAM   Compute the element forces for a beam element in the LCS.
%
%   Forces = FORCESLCS_BEAM(KeLCS,UeLCS,DLoadLCS,L,TLoadLCS,A,E,alpha,Iyy,Izz,hy,hz)
%   Forces = FORCESLCS_BEAM(KeLCS,UeLCS,DLoadLCS,L)
%   Forces = FORCESLCS_BEAM(KeLCS,UeLCS)
%       computes the element forces for the beam element in the local 
%       coordinate system (algebraic convention).
%
%   [Forces,dForcesdx] = FORCESLCS_BEAM(KeLCS,UeLCS,DLoadLCS,L,[],[],[],...
%                                 [],[],[],[],[],dKeLCSdx,dUeLCSdx,dDLoadLCSdx,dLdx)
%   [Forces,dForcesdx] = FORCESLCS_BEAM(KeLCS,UeLCS,[],[],[],[],[],[],[],...
%                                       [],[],[],dKeLCSdx,dUeLCSdx,dDLoadLCSdx,dLdx)
%       additionally computes the derivatives of the element forces with
%       respect to the design variables x.
%
%   KeLCS           Element stiffness matrix (12 * 12)
%   UeLCS           Displacements            (12 * nLC)
%   DLoadLCS        Distributed loads        [n1localX; n1localY; n1localZ; ...]
%   L               Beam length
%   DLoadLCS        Distributed loads       [n1localX; n1localY; n1localZ; ...]
%                                                                   (6 * 1)
%   TLoadLCS        Temperature load        [dTm; dTy; dTz] (3 * 1)
%   A               Cross-sectional area    
%   E               Young's modulus
%   alpha           Linear thermal expansion coefficient
%   Iyy             Area moment of inertia for bending around local y-axis
%   Izz             Area moment of inertia for bending around local z-axis
%   hy              Profile height (in local y-direction)
%   hz              Profile height (in local z-direction)
%   dKeLCSdx        Element stiffness matrix derivatives    (CELL(nVar,1))
%   dUeLCSdx        Displacements derivatives               (SIZE(UeLCS) * nVar)
%   dDLoadLCSdx     Distributed loads derivatives           (SIZE(DLoadLCS) * nVar)
%   dLdx            Beam length derivatives                 (1 * nVar)
%   Forces          Element forces                  [N; Vy; Vz; T; My; Mz] (12 * nLC)
%   dForcesdx       Element forces derivatives      (12 * nLC * nVar)
%
%   See also FORCES_BEAM, FORCESLCS_TRUSS, ELEMFORCES

% David Dooms, Wouter Dillen
% October 2008, April 2017

if nargin<3, DLoadLCS = []; end
if nargin<4, L = []; end
if nargin<5, TLoadLCS = []; end
if nargin<6, A = []; end
if nargin<7, E = []; end
if nargin<8, alpha = []; end
if nargin<9, Iyy = []; end
if nargin<10, Izz = []; end
if nargin<11, hy = []; end
if nargin<12, hz = []; end
if nargin<13, dKeLCSdx = []; end
if nargin<14, dUeLCSdx = []; end
if nargin<15, dDLoadLCSdx = []; end
if nargin<16, dLdx = []; end

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
    
% subtract equivalent nodal forces (distributed loads) from element forces
if ~isempty(DLoadLCS)
    nDLoadsOnElement = size(DLoadLCS,3);
    for i=1:nDLoadsOnElement
        [eqForces_i,deqForcesdx_i] = loadslcs_beam(DLoadLCS(:,:,i),L,dDLoadLCSdx(:,:,:,i),dLdx);
        Forces = Forces - eqForces_i;
        dForcesdx = dForcesdx - deqForcesdx_i;
    end
end

if ~isempty(TLoadLCS)
    dTm=TLoadLCS(1,:);          % Average element temperature
    dTy=TLoadLCS(2,:);          % Temperature gradient in local y direction
    dTz=TLoadLCS(3,:);          % Temperature gradient in local z direction

    Flcs=zeros(12,size(dTm,2));
    % Axial forces
    Flcs(1,:)=-E*A*alpha*dTm;
    Flcs(7,:)= E*A*alpha*dTm;
    % Bending moments Mz
    if any(TLoadLCS(2,:))
        Flcs(6,:)=  E*Izz*alpha*dTy/hy;
        Flcs(12,:)=-E*Izz*alpha*dTy/hy;
    end
    % Bending moments My
    if any(TLoadLCS(3,:))
        Flcs(5,:)= -E*Iyy*alpha*dTz/hz;
        Flcs(11,:)= E*Iyy*alpha*dTz/hz;
    end
    Forces=Forces-Flcs;
end
            
