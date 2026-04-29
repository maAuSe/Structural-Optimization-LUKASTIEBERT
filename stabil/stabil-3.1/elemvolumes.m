function [V,dVdx] = elemvolumes(Nodes,Elements,Types,Sections,dNodesdx,dSectionsdx)

%ELEMVOLUMES   Compute element volumes.
%
%   V = ELEMVOLUMES(Nodes,Elements,Types,Sections)
%   computes the volume for all elements.
%   
%   [V,dVdx] = ELEMVOLUMES(Nodes,Elements,Types,Sections,dNodesdx,dSectionsdx)
%   computes the volume for all elements, as well as the derivatives of the
%   volume with respect to the design variables x.
%
%   Nodes       Node definitions          [NodID x y z]
%   Elements    Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Types       Element type definitions  {TypID EltName Option1 ... }
%   Sections    Section definitions       [SecID SecProp1 SecProp2 ...]
%   dNodesdx    Node definitions derivatives     (SIZE(Nodes) * nVar)
%   dSectionsdx Section dervinitions derivatives (SIZE(Sections) * nVar)
%   V           Element volumes              (nElem * 1)
%   dVdx        Element volumes derivatives  (nElem * nVar)
%
%   See also VOLUME_BEAM, VOLUME_TRUSS, ELEMSIZES.

% Mattias Schevenels, Wouter Dillen
% March 2017, December 2017

% preprocessing
if nargin<4, Sections = []; end
if nargin<5, dNodesdx = []; end
if nargin<6, dSectionsdx = []; end

nVar = 0;
if nargout>1 && (~isempty(dNodesdx) || ~isempty(dSectionsdx))
    nVar = max(size(dNodesdx,3),size(dSectionsdx,3));
end

if nVar==0 || isempty(dNodesdx), dNodesdx = zeros([size(Nodes),nVar]); end
if nVar==0 || isempty(dSectionsdx), dSectionsdx = zeros([size(Sections),nVar]); end


V = zeros(size(Elements,1),1);
dVdx = zeros([length(V),nVar]);

for iType=1:size(Types,1)
    loc2=find(Elements(:,2)==cell2mat(Types(iType,1)));
    
    if ~isempty(loc2)
        Type=Types{iType,2};
        
        for i=1:length(loc2)
            iElem = loc2(i);
            
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
            
            % Section
            SecID=Elements(iElem,3);
            loc=find(Sections(:,1)==SecID);
            if isempty(loc)
                error('Section %i is not defined.',SecID)
            elseif length(loc)>1
                error('Section %i is multiply defined.',SecID)
            end
            Section=Sections(loc,2:end);
            dSectiondx=dSectionsdx(loc,2:end,:);
            
            [V(iElem),dVdx(iElem,:)] = eval(['volume_' Type '(Node,Section,dNodedx,dSectiondx)']);
        end
    end
end

