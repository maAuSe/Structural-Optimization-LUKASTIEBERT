function [ForcesLCS,ForcesGCS,dForcesLCSdx,dForcesGCSdx]=forces_truss(Node,Section,Material,UeGCS,DLoad,TLoad,Options,dNodedx,dSectiondx,dUeGCSdx,dDLoaddx)

%FORCES_TRUSS   Compute the element forces for a truss element.
%
%   [ForcesLCS,ForcesGCS] = FORCES_TRUSS(Node,Section,Material,UeGCS,[],TLoad)
%   [ForcesLCS,ForcesGCS] = FORCES_TRUSS(Node,Section,Material,UeGCS)
%       computes the element forces for the truss element in the local and the  
%       global coordinate system (algebraic convention).
%
%   [ForcesLCS,ForcesGCS,dForcesLCSdx,dForcesGCSdx]
%           = FORCES_TRUSS(Node,Section,Material,UeGCS,[],TLoad,[],dNodedx,...
%                                                                dSectiondx,dUeGCSdx)
%       additionally computes the derivatives of the element forces with
%       respect to the design variables x.
%
%   Node            Node definitions           [x y z] (2 * 3)
%   Section         Section definition         [A]
%   Material        Material definition        [E]
%   UeGCS           Displacements              (6 * nLC)
%   TLoad           Temperature load           [dTm]
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
%   See also FORCESLCS_TRUSS, ELEMFORCES.

% David Dooms, Wouter Dillen
% October 2008, April 2017

if nargin<6, TLoad = []; end
if nargin<8, dNodedx = []; end
if nargin<9, dSectiondx = []; end
if nargin<10, dUeGCSdx = []; end

nVar = 0;
if nargout>2 && (~isempty(dNodedx) || ~isempty(dSectiondx) || ~isempty(dUeGCSdx))
    nVar = max([size(dNodedx,3),size(dSectiondx,3),size(dUeGCSdx,3)]);
end 

if nVar==0 || isempty(dNodedx), dNodedx = zeros([size(Node),nVar]); end
if nVar==0 || isempty(dSectiondx), dSectiondx = zeros([size(Section),nVar]); end
if nVar==0 || isempty(dUeGCSdx), dUeGCSdx = zeros(size(UeGCS,1),size(UeGCS,2),nVar); end

dForcesLCSdx=zeros(12,size(dUeGCSdx,2),nVar);
dForcesGCSdx=zeros(12,size(dUeGCSdx,2),nVar);


% Element length
L=norm(Node(2,:)-Node(1,:)); 
dLdx=(permute(dNodedx(2,:,:)-dNodedx(1,:,:),[3 2 1])*(Node(2,:)-Node(1,:)).')/L;

% Material
E=Material(1);

% Section
A=Section(1);
dAdx=dSectiondx(:,1,:);
alpha=[];
if length(Material)>3
    alpha=Material(4); % Material=[E nu rho alpha]
elseif length(Material)>2
    alpha=Material(3); % Material=[E nu alpha]
end
    
[t,dtdx]=trans_truss(Node,dNodedx);
T=blkdiag(t,t);
dTdx = [dtdx  zeros(size(dtdx));
        zeros(size(dtdx)) dtdx];
        
% UeLCS
UeLCS=T*UeGCS;
dUeLCSdx=zeros(size(dUeGCSdx));
for n=1:nVar
    dUeLCSdx(:,:,n)=dTdx(:,:,n)*UeGCS + T*dUeGCSdx(:,:,n);
end

% TLoadLCS
if isempty(TLoad)
    dTm=zeros(1,size(UeGCS,2));
else
    dTm=TLoad(1,:);
end
    
% KeLCS
[KeLCS,~,dKeLCSdx]=kelcs_truss(L,A,E,[],{},dLdx,dAdx);
    
% ForcesLCS
if isempty(alpha)
    [forcesLCS1,dforcesLCS1dx]=forceslcs_truss(KeLCS,UeLCS,[],[],[],[],dKeLCSdx,dUeLCSdx);
else
    [forcesLCS1,dforcesLCS1dx]=forceslcs_truss(KeLCS,UeLCS,dTm,A,E,alpha,dKeLCSdx,dUeLCSdx);
end

% ForcesGCS
forcesGCS1=T.'*forcesLCS1;
dforcesGCS1dx=zeros(size(dforcesLCS1dx));
for n=1:nVar
    dforcesGCS1dx(:,:,n)=dTdx(:,:,n).'*forcesLCS1+T.'*dforcesLCS1dx(:,:,n);
end
    
% reorder
ForcesLCS=zeros(12,size(UeGCS,2));
ForcesLCS([1,7],:)=forcesLCS1([1,4],:);
ForcesGCS=zeros(12,size(UeGCS,2));
ForcesGCS([1:3,7:9],:)=forcesGCS1;
if nVar>0
    dForcesLCSdx([1,7],:,:)=dforcesLCS1dx([1,4],:,:);
    dForcesGCSdx([1:3,7:9],:,:)=dforcesGCS1dx;
end
