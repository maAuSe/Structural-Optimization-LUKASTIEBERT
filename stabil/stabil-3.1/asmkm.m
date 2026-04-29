function [K,M,dKdx] = asmkm(Nodes,Elements,Types,Sections,Materials,DOF,dNodesdx,dSectionsdx)

%ASMKM   Assemble stiffness and mass matrix.
%
%   [K,M] = ASMKM(Nodes,Elements,Types,Sections,Materials,DOF)
%    K    = ASMKM(Nodes,Elements,Types,Sections,Materials,DOF)
%    K    = ASMKM(Nodes,Elements,Types,Sections,Materials)
%   assembles the stiffness and the mass matrix using the finite element method.
%
%   [K,~,dKdx] = ASMKM(Nodes,Elements,Types,Sections,Materials,DOF,dNodesdx,dSectionsdx)
%   assembles the stiffness matrix using the finite element method and
%   additionally computes the derivatives of the stiffness matrix with
%   respect to the design variables x. The derivatives of the mass matrix
%   have not yet been implemented.
%
%   Nodes       Node definitions          [NodID x y z]
%   Elements    Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Types       Element type definitions  {TypID EltName Option1 ...}
%   Sections    Section definitions       [SecID SecProp1 SecProp2 ...]
%   Materials   Material definitions      [MatID MatProp1 MatProp2 ...]
%   DOF         Degrees of freedom              (nDOF * 1)
%   dNodesdx    Node definitions derivatives    (SIZE(Nodes) * nVar)
%   dSectionsdx Section definitions derivatives (SIZE(Sections) * nVar)
%   K           Stiffness matrix                (nDOF * nDOF)
%   M           Mass matrix                     (nDOF * nDOF)
%   dKdx        Stiffness matrix derivatives    (CELL(nVar,1))
%
%   See also KE_TRUSS, KE_BEAM.

if nargin<7, dNodesdx = []; end
if nargin<8, dSectionsdx = []; end

nVar = 0;
if nargout>2 && (~isempty(dNodesdx) || ~isempty(dSectionsdx))
    nVar = max(size(dNodesdx,3),size(dSectionsdx,3));
end
  
if ~isempty(dNodesdx) && ~isempty(dSectionsdx) && (size(dNodesdx,3)~=size(dSectionsdx,3))
    error('SIZE(dNodesdx,3) and SIZE(dSectionsdx,3) must be identical.');
end

if nVar==0 || isempty(dNodesdx), dNodesdx = zeros([size(Nodes),nVar]); end
if nVar==0 || isempty(dSectionsdx), dSectionsdx = zeros([size(Sections),nVar]); end

if nargout>3
    error('Sensitivities for the mass matrix have not been implemented yet.')
end


% initialize DOF
[aDOF,elemInd] = getdof(Elements,Types);
if nargin < 6 || isempty(DOF)
    DOF = aDOF;
    nDOF=length(DOF);
    adofInd = (1:nDOF).';
else 
    DOF = DOF(:);
    nDOF=length(DOF);
    [~,adofInd] = ismember(round(aDOF*100),round(DOF*100));
end

nElem=size(Elements,1);

%% Perform checks on input matrices
checkunique(Nodes(:,1),'Node');
checkunique(cell2mat(Types(:,1)),'Element type');
checkunique(DOF,'Degree of freedom');
checkunique(Sections(:,1),'Section');
checkunique(Materials(:,1),'Material');

%% Element stiffness matrix functions
kefunc = cellfun(@(x) ['ke_' x],Types(:,2),'UniformOutput',0);
kefunc = cellfun(@str2func,kefunc,'UniformOutput',0);

%% Link Elements with Types, Sections, Materials
Elements(:,2)=linkandcheck(round(Elements(:,2)),round(cell2mat(Types(:,1))),'Element type');
Elements(:,3)=linkandcheck(round(Elements(:,3)),round(Sections(:,1)),'Section');
Elements(:,4)=linkandcheck(round(Elements(:,4)),round(Materials(:,1)),'Material');

[~,Elements(:,5:end)] = ismember(round(Elements(:,5:end)),round(Nodes(:,1)));


%% Assemble in blocks to save memory
nElemMax=2000; % number of elements per block
nBlock=ceil(nElem/nElemMax);
Blocks = [(0:nElemMax:nElemMax*(nBlock-1)),nElem];

%% Allocate output
K=sparse(nDOF,nDOF);
if nargout>1, M=sparse(nDOF,nDOF); end
if nargout>2
    dKdx=cell(nVar,1);
    for n=1:nVar, dKdx{n}=sparse(nDOF,nDOF); end
end

%% Main loop over elements per block
for iBlock=1:nBlock
Blstart = Blocks(iBlock)+1;
Blend = Blocks(iBlock+1);
for iElem = Blstart:Blend
    
    if size(Types,2)<3
        Options={};
    else
        Options=Types{Elements(iElem,2),3};
    end
    
    Section=Sections(Elements(iElem,3),2:end);
    dSectiondx=dSectionsdx(Elements(iElem,3),2:end,:);
        
    Material=Materials(Elements(iElem,4),2:end); 
    
    NodeNum=Elements(iElem,5:end);
    Node=NaN(length(NodeNum),3);
    Node(NodeNum~=0,:) = Nodes(NodeNum(NodeNum~=0),2:end);
    dNodedx=NaN([size(Node),nVar]);
    dNodedx(NodeNum~=0,:,:) = dNodesdx(NodeNum(NodeNum~=0),2:end,:);
     
    inddofelem = adofInd(elemInd{iElem}).';    
    utilind=find(inddofelem);
    nutilind=length(utilind);
    inddofelem = inddofelem(ones(1,nutilind),utilind);
    indi{iElem}=reshape(inddofelem.',1,nutilind^2);
    indj{iElem}=reshape(inddofelem,1,nutilind^2);
    
    try        
    if nargout==2         % stiffness and mass
        [Ke,Me]=kefunc{Elements(iElem,2)}(Node,Section,Material,Options);
        sM{iElem}=reshape(Me(utilind,utilind),1,nutilind^2);
        sK{iElem}=reshape(Ke(utilind,utilind),1,nutilind^2);
    else                                % only stiffness
        [Ke,~,dKedx]=kefunc{Elements(iElem,2)}(Node,Section,Material,Options,dNodedx,dSectiondx);
        sK{iElem}=reshape(Ke(utilind,utilind),1,nutilind^2);
        for n=1:nVar
            dsKdx{(iElem-1)*nVar+n}=reshape(dKedx{n}(utilind,utilind),1,nutilind^2);
        end
    end
    catch lerr
        error('Element on line %i: %s',iElem,lerr.message)
    end
    
    
end
indi=cell2mat(indi);
indj=cell2mat(indj);
sK=cell2mat(sK);
K0=sparse(indi,indj,sK,nDOF,nDOF); clear sK;
K=K+K0;
for n=1:nVar
    tmp = cell2mat(dsKdx(n:nVar:end));
    if any(tmp~=0)
        dK0dx=sparse(indi,indj,tmp,nDOF,nDOF);
        dKdx{n}=dKdx{n}+dK0dx;
    end
end
if nVar>0, clear dsKdx; end

if nargout==2  % stiffness and mass
    sM=cell2mat(sM);
    M0=sparse(indi,indj,sM,nDOF,nDOF); clear sM;
	M=M+M0;
end

clear indi indj;
end

%% Make the matrix exactly symmetric
K=(K+K.')/2;
for n=1:nVar
    dKdx{n}=(dKdx{n}+dKdx{n}.')/2;
end
if nargout==2
    M=(M+M.')/2; 
end

end