function [Sn,Sn2]=nodalstress(Nodes,Elements,Types,Se)

%NODALSTRESS   Compute the nodal stresses from the element solution.
%
%   [Sn,Sn2] = nodalstress(Nodes,Elements,Types,Se)
%   computes the nodal stresses from the element solution.
%
%   Nodes      Node definitions          [NodID x y z]
%   Elements   Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   Sections   Section definitions       [SecID SecProp1 SecProp2 ...]
%   Se         Element stresses in corner nodes IJKL and 
%              at top/mid/bot of shell (nElem * 72)   
%              72 = 6 stresscomp. * 4 nodes * 3 locations
%                                        [sxx syy szz sxy syz sxz ...]
%   Sn         Nodal stress in corner nodes IJKL and 
%              at top/mid/bot of shell (nElem * 72)   
%              72 = 6 stresscomp. * 4 nodes * 3 locations
%                                        [sxx syy szz sxy syz sxz ...]
%   Sn2        Nodal stresses per node 
%              at top/mid/bot of shell (nNodes * 19)   
%              (19 = NodID + 3 locations * 6 scomp.)
%                                        [NodID sxx syy szz sxy syz sxz ...]
%   See also ELEMSTRESS.

% Miche Jansen
% 2010
nNode = size(Nodes,1);
Sn = zeros(size(Se));
if nargout > 1
Sn2 = zeros(nNode,19,size(Se,3));
end

for iNode=1:nNode
    Stress = [];
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
                Stress = [Stress;Se(i1,[(1:6) (25:30) (49:54)]+6*(i2-1),:)];
                inds = [inds;i1,i2];
                case {'shell6'}
                if i2 <= 3
                Stress = [Stress;Se(i1,[(1:6) (25:30) (49:54)]+6*(i2-1),:)];
                inds = [inds;i1,i2];
                elseif i2 <=6
                Stress = [Stress;(Se(i1,[(1:6) (25:30) (49:54)]+6*(i2-4),:)+Se(i1,[(1:6) (25:30) (49:54)]+6*(i2-3),:))/2];
                elseif i2==6
                Stress = [Stress;(Se(i1,[(1:6) (25:30) (49:54)]+12,:)+Se(i1,[(1:6) (25:30) (49:54)],:))/2];   
                end                
                case {'shell8'}
                if i2 <= 4
                Stress = [Stress;Se(i1,[(1:6) (25:30) (49:54)]+6*(i2-1),:)];
                inds = [inds;i1,i2];
                elseif i2 <=7
                Stress = [Stress;(Se(i1,[(1:6) (25:30) (49:54)]+6*(i2-5),:)+Se(i1,[(1:6) (25:30) (49:54)]+6*(i2-4),:))/2];
                elseif i2==8
                Stress = [Stress;(Se(i1,[(1:6) (25:30) (49:54)]+18,:)+Se(i1,[(1:6) (25:30) (49:54)],:))/2];   
                end
                case {'solid8'}
                Stress = [Stress;repmat(Se(i1,(1:6)+6*(i2-1),:),1,3)];
                inds = [inds;i1,i2];   
                case {'solid20'}
                if i2 <= 8
                Stress = [Stress;repmat(Se(i1,(1:6)+6*(i2-1),:),1,3)];
                inds = [inds;i1,i2];
                end
            end
        end
    end
    
    if ~isempty(Stress)
        Stress = mean(Stress,1);
        for iElem=1:size(inds,1)
            Sn(inds(iElem,1),6*(inds(iElem,2)-1)+[(1:6) (25:30) (49:54)],:) = Stress;
        end   
        if nargout > 1
        Sn2(iNode,1,:) = Node;
        Sn2(iNode,2:end,:) = Stress;
        end
    end
end
%Sn2 = Sn2(~isnan(Sn2));
end