function [FeLCS]=elemshellf(Elements,Sections,SeLCS)

%ELEMSHELLF   Compute the shell forces and moments per unit length for shell elements
%             (element solution).
%
%   [FeLCS] = elemshellf(Elements,Sections,SeLCS)
%   computes the element forces/moments per unit length in the local coordinate system.
%
%   Elements   Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Sections   Section definitions       [SecID SecProp1 SecProp2 ...]
%   SeLCS      Element stresses in LCS in corner nodes IJKL and
%              at top/mid/bot of shell (nElem * 72)
%                                        [sxx syy szz sxy syz sxz]
%   FeLCS      Element forces/moments per unit length (nElem * 32)
%              (32 = 8 forces * 4 corner nodes)
%                                        [Nx Ny Nxy Mx My Mxy Vx Vy]
%   See also ELEMSTRESS.

% Miche Jansen
% 2010


%
nElem=size(Elements,1);
FeLCS = zeros(nElem,32,size(SeLCS,3));

for iElem=1:nElem

    SecID=Elements(iElem,3);
    loc=find(Sections(:,1)==SecID);
    if isempty(loc)
        error('Section %i is not defined.',SecID)
    elseif length(loc)>1
        error('Section %i is multiply defined.',SecID)
    end
    Section=Sections(loc,2:end);

    h=Section(1); % werkt enkel voor constante dikte

    for iNode=1:4 % formules uit Ansys theory p2-14
    FeLCS(iElem,(1:8)+8*(iNode-1),:) = h/6*[SeLCS(iElem,1+6*(iNode-1),:)+4*SeLCS(iElem,25+6*(iNode-1),:)+SeLCS(iElem,49+6*(iNode-1),:) ...
                                       SeLCS(iElem,2+6*(iNode-1),:)+4*SeLCS(iElem,26+6*(iNode-1),:)+SeLCS(iElem,50+6*(iNode-1),:) ...
                                       SeLCS(iElem,4+6*(iNode-1),:)+4*SeLCS(iElem,28+6*(iNode-1),:)+SeLCS(iElem,52+6*(iNode-1),:) ...
                                       h*(SeLCS(iElem,1+6*(iNode-1),:)-SeLCS(iElem,49+6*(iNode-1),:))/2 ...
                                       h*(SeLCS(iElem,2+6*(iNode-1),:)-SeLCS(iElem,50+6*(iNode-1),:))/2 ...
                                       h*(SeLCS(iElem,4+6*(iNode-1),:)-SeLCS(iElem,52+6*(iNode-1),:))/2 ...
                                       SeLCS(iElem,6+6*(iNode-1),:)+4*SeLCS(iElem,30+6*(iNode-1),:)+SeLCS(iElem,54+6*(iNode-1),:) ...
                                       SeLCS(iElem,5+6*(iNode-1),:)+4*SeLCS(iElem,29+6*(iNode-1),:)+SeLCS(iElem,53+6*(iNode-1),:)];
    end
end
end