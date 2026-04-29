function [FnLCS,FnLCS2]=nodalshellf(Nodes,Elements,Types,FeLCS)

%NODALSHELLF   Compute the nodal shell forces/moments per unit length 
%                                                    from the element solution.
%
%   [FnLCS,FnLCS2] = nodalshellf(Nodes,Elements,Types,FeLCS) computes the 
%   nodal forces from the element solution.
%
%   Nodes      Node definitions          [NodID x y z]
%   Elements   Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   Sections   Section definitions       [SecID SecProp1 SecProp2 ...]
%   FeLCS      Element forces/moments per unit length in corner nodes IJKL 
%              (nElem * 32)              [Nx Ny Nxy Mx My Mxy Vx Vy]
%   FnLCS      Nodal forces/moments per unit length in corner nodes IJKL 
%              (nElem * 32)              [Nx Ny Nxy Mx My Mxy Vx Vy]
%   FnLCS2      Nodal forces per node (nNodes * 9)   
%              (9 = NodID + 8 fcomp.)
%                                        [NodID Nx Ny Nxy Mx My Mxy Vx Vy]
%   See also ELEMSHELLF.

% Miche Jansen
% 2010
nNode = size(Nodes,1);
FnLCS = zeros(size(FeLCS));
if nargout > 1
FnLCS2 = zeros(nNode,9,size(FeLCS,3));
end

for iNode=1:nNode
    shellf = [];
    inds = [];
    Node = Nodes(iNode,1);
    [indi,indj]=find(Elements(:,5:end)==Node);
    if ~isempty(indi)
        teller = 0;
        for iElem=1:length(indi)
            i1 =indi(iElem); 
            i2 =indj(iElem);
            
            TypID=Elements(i1,2);
            loc=find(cell2mat(Types(:,1))==TypID);
            if isempty(loc)
            error('Element type %i is not defined.',TypID)
            elseif length(loc)>1
            error('Element type %i is multiply defined.',TypID)
            end
    
            Type=Types{loc,2};
    
            switch lower(Type)
                case {'beam'}
                
                case {'truss'}
                    
                case {'shell4'}
                shellf = [shellf;FeLCS(i1,(1:8)+8*(i2-1),:)];
                inds = [inds;i1,i2];
                case {'shell6'}
                if i2 <= 3
                shellf = [shellf;FeLCS(i1,(1:8)+8*(i2-1),:)];
                inds = [inds;i1,i2];
                elseif i2 <=6
                shellf = [shellf;(FeLCS(i1,(1:8)+8*(i2-4),:)+FeLCS(i1,(1:8)+8*(i2-3),:))/2];
                elseif i2==6
                shellf = [shellf;(FeLCS(i1,(1:8)+16,:)+FeLCS(i1,(1:8),:))/2];   
                end
                case {'shell8'}
                if i2 <= 4
                shellf = [shellf;FeLCS(i1,(1:8)+8*(i2-1),:)];
                inds = [inds;i1,i2];
                elseif i2 <=7
                shellf = [shellf;(FeLCS(i1,(1:8)+8*(i2-5),:)+FeLCS(i1,(1:8)+8*(i2-4),:))/2];
                elseif i2==8
                shellf = [shellf;(FeLCS(i1,(1:8)+24,:)+FeLCS(i1,(1:8),:))/2];   
                end
            end
        end
    end
    
    if ~isempty(shellf)
        shellf = mean(shellf,1);
        for iElem=1:size(inds,1)
            FnLCS(inds(iElem,1),(1:8)+8*(inds(iElem,2)-1),:) = shellf;
        end 
        if nargout > 1
        FnLCS2(iNode,1,:) = Node;
        FnLCS2(iNode,2:end,:) = shellf;
        end
    end
end
%FnLCS = FnLCS(~isnan(FnLCS));
end