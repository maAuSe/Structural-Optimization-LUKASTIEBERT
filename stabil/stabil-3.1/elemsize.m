function [S,dSdx] = elemsize(Nodes,Elements,Types,dNodesdx)

%ELEMSIZE   Compute element length/area/volume.
%
%   S = ELEMSIZE(Nodes,Elements,Types)
%   computes the size of all elements, depending on element type. For a 1D 
%   line element it computes length, for a 2D plate element it computes
%   area, and for a 3D solid element it computes volume.
%   
%   [S,dSdx] = ELEMSIZE(Nodes,Elements,Types,dNodesdx)
%   additionaly computes the derivatives of the size with respect to the 
%   design variables x.
%
%   Nodes       Node definitions                [NodID x y z]
%   Elements    Element definitions             [EltID TypID SecID MatID n1 n2 ...]
%   Types       Element type definitions        {TypID EltName Option1 ... }
%   dNodesdx    Node definitions derivatives 	(SIZE(Nodes) * nVar)
%   S           Element sizes
%   dSdx        Element sizes derivatives
%
%   See also SIZE_BEAM, SIZE_TRUSS, ELEMVOLUMES.

% Wouter Dillen
% December 2017

if nargin<4, dNodesdx = []; end

nVar = 0;
if nargout>1 && ~isempty(dNodesdx)
    nVar = size(dNodesdx,3);
end

if nVar==0 || isempty(dNodesdx), dNodesdx = zeros([size(Nodes),nVar]); end


S = zeros(size(Elements,1),1);
dSdx = zeros([length(S),nVar]);

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
            
            if nargout>1
                [S(iElem),dSdx(iElem,:)] = eval(['size_' Type '(Node,dNodedx)']);
            else
                S(iElem) = eval(['size_' Type '(Node)']);
            end
            
        end
    end
end

