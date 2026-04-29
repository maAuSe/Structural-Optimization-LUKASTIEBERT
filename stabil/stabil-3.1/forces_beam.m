function [ForcesLCS,ForcesGCS,dForcesLCSdx,dForcesGCSdx]=forces_beam(Node,Section,Material,UeGCS,DLoad,TLoad,Options,dNodedx,dSectiondx,dUeGCSdx,dDLoaddx)

%FORCES_BEAM   Compute the element forces for a beam element.
%
%   [ForcesLCS,ForcesGCS] = FORCES_BEAM(Node,Section,Material,UeGCS,DLoad,TLoad,Options)
%   [ForcesLCS,ForcesGCS] = FORCES_BEAM(Node,Section,Material,UeGCS,DLoad)
%   [ForcesLCS,ForcesGCS] = FORCES_BEAM(Node,Section,Material,UeGCS)
%       computes the element forces for the beam element in the local and 
%       the global coordinate system (algebraic convention).
%
%   [ForcesLCS,ForcesGCS,dForcesLCSdx,dForcesGCSdx]
%           = FORCES_BEAM(Node,Section,Material,UeGCS,DLoad,TLoad,Options,dNodedx,...
%                                                          dSectiondx,dUeGCSdx,dDLoaddx)
%           = FORCES_BEAM(Node,Section,Material,UeGCS,[],[],Options,dNodedx,...
%                                                                   dSectiondx,dUeGCSdx)
%       additionally computes the derivatives of the element forces with
%       respect to the design variables x.
%
%   Node            Node definitions           [x y z] (3 * 3)
%   Section         Section definition         [A ky kz Ixx Iyy Izz]
%   Material        Material definition        [E nu]
%   UeGCS           Displacements              (12 * nLC)
%   DLoad           Distributed loads       [n1globalX; n1globalY; n1globalZ; ...] 
%                                                                        (6 * 1)
%   TLoad           TLoad                      [dTm; dTy; dTz] (3 * 1)
%   Options         Element options            {Option1 Option2 ...}
%   dNodedx         Node definitions derivatives        (SIZE(Node) * nVar)
%   dSectiondx      Section definitions derivatives     (SIZE(Section) * nVar)
%   dUeGCSdx        Displacements derivatives           (SIZE(UeGCS) * nVar)
%   dDLoaddx        Distributed loads derivatives       (SIZE(DLoad) * nVar)
%   ForcesLCS       Element forces in the LCS           (12 * nLC)
%   ForcesGCS       Element forces in the GCS           (12 * nLC)
%   dForcesLCSdx    Element forces derivatives in LCS   (12 * nLC * nVar)
%   dForcesGCSdx    Element forces derivatives in GCS   (12 * nLC * nVar)
%
%   See also FORCESLCS_BEAM, ELEMFORCES.

% David Dooms, Wouter Dillen
% October 2008, April 2017

if nargin<5, DLoad = []; end
if nargin<6, TLoad = []; end
if nargin<8, dNodedx = []; end
if nargin<9, dSectiondx = []; end
if nargin<10, dUeGCSdx = []; end
if nargin<11, dDLoaddx = []; end

nVar = 0;
if nargout>2 && (~isempty(dNodedx) || ~isempty(dSectiondx) || ~isempty(dUeGCSdx) || ~isempty(dDLoaddx))
    nVar = max([size(dNodedx,3),size(dSectiondx,3),size(dUeGCSdx,3),size(dDLoaddx,3)]);
end 

if nVar==0 || isempty(dNodedx), dNodedx = zeros([size(Node),nVar]); end
if nVar==0 || isempty(dSectiondx), dSectiondx = zeros([size(Section),nVar]); end
if nVar==0 || isempty(dUeGCSdx), dUeGCSdx = zeros(size(UeGCS,1),size(UeGCS,2),nVar); end
if nVar==0 || isempty(dDLoaddx), dDLoaddx = zeros(size(DLoad,1),size(DLoad,2),nVar,size(DLoad,3)); end

%dForcesLCSdx=zeros(size(dUeGCSdx));
dForcesGCSdx=zeros(size(dUeGCSdx));


% Element length
L=norm(Node(2,:)-Node(1,:)); 
dLdx=(permute(dNodedx(2,:,:)-dNodedx(1,:,:),[3 2 1])*(Node(2,:)-Node(1,:)).')/L;

% Material
E=Material(1);
nu=Material(2);
alpha=[];
if length(Material)>3
    alpha=Material(4); % Material=[E nu rho alpha]
elseif length(Material)>2
    alpha=Material(3); % Material=[E nu alpha]
end

% Section
A=Section(1);
ky=Section(2);
kz=Section(3);
Ixx=Section(4);
Iyy=Section(5);
Izz=Section(6);
dAdx=dSectiondx(:,1,:);
dkydx=dSectiondx(:,2,:);
dkzdx=dSectiondx(:,3,:);
dIxxdx=dSectiondx(:,4,:);
dIyydx=dSectiondx(:,5,:);
dIzzdx=dSectiondx(:,6,:);

% Transformation matrix
[t,dtdx] = trans_beam(Node,dNodedx);
T=blkdiag(t,t,t,t);
dTdx = [dtdx                zeros(size(dtdx)) 	zeros(size(dtdx))   zeros(size(dtdx));
        zeros(size(dtdx))   dtdx                zeros(size(dtdx))   zeros(size(dtdx));
        zeros(size(dtdx))   zeros(size(dtdx))   dtdx                zeros(size(dtdx));
        zeros(size(dtdx))   zeros(size(dtdx))   zeros(size(dtdx))   dtdx];

% UeLCS
UeLCS=T*UeGCS;
dUeLCSdx=zeros(size(dUeGCSdx));
for n=1:nVar
    dUeLCSdx(:,:,n)=dTdx(:,:,n)*UeGCS + T*dUeGCSdx(:,:,n);
end
    
% DLoadLCS
if isempty(DLoad)
    DLoadLCS=zeros(6,size(UeGCS,2));
    dDLoadLCSdx=zeros(size(DLoadLCS,1),size(DLoadLCS,2),nVar);
else
    T=blkdiag(t,t);
    dTdx = [dtdx  zeros(size(dtdx));
            zeros(size(dtdx)) dtdx];
    [DLoadLCS,dDLoadLCSdx]=dloadgcs2lcs(T,DLoad,dTdx,dDLoaddx);
end

% TLoadLCS
if isempty(TLoad)
    TLoadLCS=zeros(3,size(UeGCS,2));
else
    TLoadLCS=TLoad;
end
if any(TLoadLCS(2,:))
    try
        hy=Section(7)+Section(8);
    catch
        error('Properly define section dimensions yt and yb when applying a temperature gradient.');
    end
else
    hy=[];
end
if any(TLoadLCS(3,:))
    try
        hz=Section(9)+Section(10);
    catch
        error('Properly define section dimensions zt and zb when applying a temperature gradient.');
    end
else
    hz=[];
end
    
% KeLCS
[KeLCS,~,dKeLCSdx]=kelcs_beam(L,A,ky,kz,Ixx,Iyy,Izz,E,nu,[],{},dLdx,dAdx,dkydx,dkzdx,dIxxdx,dIyydx,dIzzdx);
    
% ForcesLCS
if isempty(alpha)
    [ForcesLCS,dForcesLCSdx]=forceslcs_beam(KeLCS,UeLCS,DLoadLCS,L,[],[],[],[],[],[],[],[],dKeLCSdx,dUeLCSdx,dDLoadLCSdx,dLdx);
else
    [ForcesLCS,dForcesLCSdx]=forceslcs_beam(KeLCS,UeLCS,DLoadLCS,L,TLoadLCS,A,E,alpha,Iyy,Izz,hy,hz,dKeLCSdx,dUeLCSdx,dDLoadLCSdx,dLdx);
end

% ForcesGCS
T=blkdiag(t,t,t,t);
ForcesGCS=T.'*ForcesLCS;
dTdx = [dtdx                zeros(size(dtdx)) 	zeros(size(dtdx))   zeros(size(dtdx));
        zeros(size(dtdx))   dtdx                zeros(size(dtdx))   zeros(size(dtdx));
        zeros(size(dtdx))   zeros(size(dtdx))   dtdx                zeros(size(dtdx));
        zeros(size(dtdx))   zeros(size(dtdx))   zeros(size(dtdx))   dtdx];
for n=1:nVar
    if any(reshape(dtdx(:,:,n),[],1)~=0)
        dForcesGCSdx(:,:,n) = dTdx(:,:,n).'*ForcesLCS + T.'*dForcesLCSdx(:,:,n);
    else
        dForcesGCSdx(:,:,n) = T.'*dForcesLCSdx(:,:,n);
    end
end
    
  
