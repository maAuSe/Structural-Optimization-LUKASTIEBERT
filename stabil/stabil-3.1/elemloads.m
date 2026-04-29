function [F,dFdx] = elemloads(DLoads,Nodes,Elements,Types,DOF,dDLoadsdx,dNodesdx)

%ELEMLOADS   Equivalent nodal forces.
%
%   F = ELEMLOADS(DLoads,Nodes,Elements,Types,DOF)
%   computes the equivalent nodal forces of a distributed load 
%   (in the global coordinate system).
%
%   [F,dFdx] = ELEMLOADS(DLoads,Nodes,Elements,Types,DOF,dDLoadsdx,dNodesdx)
%   additionally computes the derivatives of the equivalent nodal forces
%   with respect to the design variables x.
%
%   DLoads     Distributed loads                [EltID n1globalX n1globalY n1globalZ ...]
%   Nodes      Node definitions                 [NodID x y z]
%   Elements   Element definitions              [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions         {TypID EltName Option1 ... }
%   DOF        Degrees of freedom               (nDOF * 1)
%   dDLoadsdx  Distributed loads derivatives    (SIZE(DLoads) * nVar)
%   dNodesdx   Node definitions derivatives     (SIZE(Nodes) * nVar)
%   F          Load vector              (nDOF * nLC)
%   dFdx       Load vector derivatives  (nDOF * nLC * nVar)
%
%   See also LOADS_TRUSS, LOADS_BEAM, NODALVALUES.

% David Dooms, Wouter Dillen
% March 2008, April 2017

if nargin<6, dDLoadsdx = []; end
if nargin<7, dNodesdx = []; end

if ~isempty(dDLoadsdx) && ~isequal(DLoads(:,1),dDLoadsdx(:,1,1,1))
    error('The first column (EltID) of DLoads and dDLoadsdx should be identical.')
end

nVar = 0;
if nargout>1 && (~isempty(dDLoadsdx) || ~isempty(dNodesdx))
    nVar = max(size(dDLoadsdx,4),size(dNodesdx,3));
end

if ~isempty(dDLoadsdx) && ~isempty(dNodesdx) && (size(dDLoadsdx,4)~=size(dNodesdx,3))
    error('SIZE(dDLoadsdx,4) and SIZE(dNodesdx,3) must be identical.');
end

if nVar==0 || isempty(dDLoadsdx), dDLoadsdx = zeros(size(DLoads,1),size(DLoads,2),size(DLoads,3),nVar); end
if nVar==0 || isempty(dNodesdx), dNodesdx = zeros([size(Nodes),nVar]); end


nDLoads=size(DLoads,1);
nterm=0;
for iDLoad=1:nDLoads
    
    % Element
    EltID=DLoads(iDLoad,1);
    loc=find(Elements(:,1)==EltID);
    if isempty(loc)
        error('Element %i is not defined.',EltID)
    elseif length(loc)>1
        error('Element %i is multiply defined.',EltID)
    end
    
    TypID=Elements(loc,2);
    NodeNum=Elements(loc,5:end);

    % Type
    loc=find(cell2mat(Types(:,1))==TypID);
    if isempty(loc)
        error('Element type %i is not defined.',TypID)
    elseif length(loc)>1
        error('Element type %i is multiply defined.',TypID)
    end
    
    Type=Types{loc,2};
    
    % Nodes
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
    
    % Fe
    DLoad=permute(DLoads(iDLoad,2:end,:),[2 3 1]);
    dDLoaddx=permute(dDLoadsdx(iDLoad,2:end,:,:),[2 3 4 1]);
    [Fe,dFedx]=eval(['loads_' Type '(DLoad,Node,dDLoaddx,dNodedx)']);
    
    PLoad(nterm+1:nterm+size(Fe,1),:)=Fe;
    dPLoaddx(nterm+1:nterm+size(Fe,1),:,:)=dFedx;      
    seldof(nterm+1:nterm+size(Fe,1),1)=eval(['dof_' Type '(NodeNum)']);
    nterm=nterm+size(Fe,1);
    
end


% F
F=nodalvalues(DOF,seldof,PLoad);
dFdx=zeros(size(F,1),size(F,2),nVar);
for n=1:nVar
    dFdx(:,:,n)=nodalvalues(DOF,seldof,dPLoaddx(:,:,n));
end

