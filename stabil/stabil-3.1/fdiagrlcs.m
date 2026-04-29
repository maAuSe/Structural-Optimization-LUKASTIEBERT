function [FdiagrLCS,dFdiagrLCSdx] = fdiagrlcs(ftype,Nodes,Elements,Types,Forces,DLoads,Points,dNodesdx,dForcesdx,dDLoadsdx)

%FDIAGRLCS   Return force diagrams in LCS.
%
%   FdiagrLCS = FDIAGRLCS(ftype,Nodes,Elements,Types,Forces,DLoads,Points)
%   FdiagrLCS = FDIAGRLCS(ftype,Nodes,Elements,Types,Forces,DLoads)
%   FdiagrLCS = FDIAGRLCS(ftype,Nodes,Elements,Types,Forces)
%       computes element force values in all interpolation points (in beam
%       convention) for beam and truss elements.
%
%   [FdiagrLCS,dFdiagrLCSdx] 
%           = FDIAGRLCS(ftype,Nodes,Elements,Types,Forces,DLoads,...
%                                           Points,dNodesdx,dForcesdx,dDLoadsdx)
%       additionally computes the derivatives of the force values with
%       respect to the design variables x.
%
%   ftype           'norm'       Normal force (in the local x-direction)
%                   'sheary'     Shear force in the local y-direction
%                   'shearz'     Shear force in the local z-direction
%                   'momx'       Torsional moment (around the local x-direction)
%                   'momy'       Bending moment around the local y-direction
%                   'momz'       Bending moment around the local z-direction
%   Nodes           Node definitions          [NodID x y z]
%   Elements        Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Types           Element type definitions  {TypID EltName Option1 ... }
%   Forces          Element forces in LCS  (beam convention) [N Vy Vz T My Mz] 
%                                                                   (nElem * 12)
%   DLoads          Distributed loads         [EltID n1globalX n1globalY n1globalZ ...]
%   Points          Points in the local coordinate system  (1 * nPoints)
%   dNodesdx        Node definitions derivatives        (SIZE(Nodes) * nVar)
%   dForcesdx       Element forces in LCS derivatives   (SIZE(Forces) * nVar)
%   dDLoadsdx       Distributed loads derivatives       (SIZE(DLoads) * nVar)
%   FdiagrLCS       Element force values at the points  (nElem * nPoints * nLC)
%   dFdiagrLCSdx    Element force values derivatives    (nElem * nPoints * nLC * nVar)
%
%   See also PLOTFORC, FDIAGRGCS_BEAM, FDIAGRGCS_TRUSS.

% Wouter Dillen
% April 2017

if nargin<6, DLoads = []; end
if nargin<7, Points = []; end
if nargin<8, dNodesdx = []; end
if nargin<9, dForcesdx = []; end
if nargin<10, dDLoadsdx = []; end

nVar = 0;
if nargout>1 && (~isempty(dNodesdx) || ~isempty(dForcesdx) || ~isempty(dDLoadsdx))
    nVar = max([size(dNodesdx,3),size(dForcesdx,4),size(dDLoadsdx,4)]);
end

if nVar==0 || isempty(dNodesdx), dNodesdx = zeros([size(Nodes),nVar]); end
if nVar==0 || isempty(dForcesdx), dForcesdx = zeros(size(Forces,1),size(Forces,2),size(Forces,3),nVar); end
if nVar==0 || isempty(dDLoadsdx), dDLoadsdx = zeros(size(DLoads,1),size(DLoads,2),size(DLoads,3),nVar); end


if ~isempty(Points)
    nPoints=length(Points);
else   
    nPoints=21;
    Points=linspace(0,1,nPoints);
end

nElem=size(Elements,1);
nLC=size(Forces,3);

FdiagrLCS=zeros(nElem,nPoints,nLC);
dFdiagrLCSdx=zeros(nElem,nPoints,nLC,nVar);


for iElem=1:nElem

    % Type
    TypID=Elements(iElem,2);
    loc=find(cell2mat(Types(:,1))==TypID);
    if isempty(loc)
        error('Element type %i is not defined.',TypID)
    elseif length(loc)>1
        error('Element type %i is multiply defined.',TypID)
    end
    
    Type=Types{loc,2};
    if not(strcmp(Type,'beam') || strcmp(Type,'truss'))
        error('Only beam and truss elements are supported.')
    end
            
    % Nodes
    NodeNum=Elements(iElem,5:end);
    Node=zeros(length(NodeNum),3);
    dNodedx=zeros(length(NodeNum),3,nVar);
    
    for iNode=1:length(NodeNum)
        loc=find(Nodes(:,1)==NodeNum(1,iNode));
        if isempty(loc)
            Node(iNode,:)=NaN;
            dNodedx(iNode,:,:)=NaN;
        elseif length(loc)>1
            error('Node %i is multiply defined.',NodeNum(1,iNode))
        else
            Node(iNode,:)=Nodes(loc,2:end);
            dNodedx(iNode,:,:)=dNodesdx(loc,2:end,:);
        end
    end
    
    % element length
    L=norm(Node(2,:)-Node(1,:));
    dLdx=(permute(dNodedx(2,:,:)-dNodedx(1,:,:),[3 2 1])*(Node(2,:)-Node(1,:)).')/L;
    
    % Forces
    Force=permute(Forces(iElem,:,:),[2 3 1]);
    dForcedx=permute(dForcesdx(iElem,:,:,:),[2 3 4 1]);

    EltID=Elements(iElem,1);

    % DLoads
    if isempty(DLoads)
        DLoads=zeros(1,7,nLC);
    end
    loc=find(DLoads(:,1)==EltID);
    if isempty(loc)
        DLoad=zeros(6,nLC);
        dDLoaddx=zeros(6,nLC,nVar);
    else
        DLoad=permute(DLoads(loc,2:end,:),[2 3 1]);
        dDLoaddx=permute(dDLoadsdx(loc,2:end,:,:),[2 3 4 1]);
    end
    
    % transform to LCS, add to FdiagrLCS
    [t,dtdx] = eval(['trans_' Type '(Node,dNodedx)']);
    T = blkdiag(t,t);
    dTdx = [dtdx  zeros(size(dtdx));
            zeros(size(dtdx)) dtdx];
    [DLoadLCS,dDLoadLCSdx] = dloadgcs2lcs(T,DLoad,dTdx,dDLoaddx);
    [FdiagrLCS(iElem,:,:),~,~,dFdiagrLCSdx(iElem,:,:,:)] = eval(['fdiagrlcs_' Type '(ftype,Force,DLoadLCS,L,Points,dForcedx,dDLoadLCSdx,dLdx)']);

end



