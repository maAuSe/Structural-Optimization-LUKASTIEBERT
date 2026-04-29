function [Spr,Vpr]=principalstress(Elements,SeGCS)

%PRINCIPALSTRESS    Compute the principal stresses and directions in shell elements.
%
%   [Spr,Vpr] = principalstress(Elements,SeGCS)
%   computes the principal stresses and directions.
%
%   Elements   Element definitions   [EltID TypID SecID MatID n1 n2...]
%   SeGCS      Element stresses in GCS in corner nodes IJKL and 
%              at top/mid/bot of shell (nElem * 72)   
%              72 = 6 stresscomp. * 4 nodes * 3 locations
%                                        [sxx syy szz sxy syz sxz ...]
%   Spr        Principal stress (nElem * 72)
%                                        [s1 s2 s3 0 0 0 ...]
%   Vpr        Principal stress directions {nElem * 12}
%
%   See also ELEMSTRESS, PLOTPRINCDIR.

nElem = size(Elements,1);
Spr = zeros(size(SeGCS));
Vpr = cell([size(SeGCS,1),12]);

for iElem=1:nElem
    
    if nnz(isnan(SeGCS(iElem,:))) == 0
        for iloc = 1:3
        for iNode = 1:4
            Se2 = SeGCS(iElem,(1:6)+6*(iNode-1)+24*(iloc-1));
            Se2 = [Se2(1),Se2(4),Se2(6);Se2(4),Se2(2),Se2(5);Se2(6),Se2(5),Se2(3)];
            [Vpr{iElem,4*(iloc-1)+iNode},Spr2]=eig(Se2);
            Spr(iElem,(1:3)+6*(iNode-1)+24*(iloc-1))=diag(Spr2);
        end
        end
    else Spr(iElem,:) = nan(1,size(SeGCS,2));
    end
end
