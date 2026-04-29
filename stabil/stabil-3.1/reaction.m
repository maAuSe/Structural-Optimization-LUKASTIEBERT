function Freac=reaction(Elements,ForcesGCS,seldof)

%REACTION   Compute the reaction forces for a degree of freedom.
%
%   Freac=reaction(Elements,ForcesGCS,seldof)
%   computes the reaction forces for a degree of freedom.
%
%   Elements   Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   ForcesGCS  Element forces in GCS (nElem * 12)
%   seldof     Selected DOF labels (kDOF * 1)
%   Freac      Reaction force 
%
%   See also ELEMFORCES.

% David Dooms
% September 2008

% Currently no wild cards allowed!!!
% Only for beams and trusses

% PREPROCESSING
seldof=seldof(:);

Freac=zeros(size(seldof));

for iDOF=1:size(seldof,1)
    NodeNum=floor(seldof(iDOF,1));
    indj=round(rem(seldof(iDOF,1),1)*100);
    loci=find(Elements(:,5)==NodeNum);
    locj=find(Elements(:,6)==NodeNum);
    Freac(iDOF,1)=Freac(iDOF,1)+sum(ForcesGCS(loci,indj));
    Freac(iDOF,1)=Freac(iDOF,1)+sum(ForcesGCS(locj,(indj+6)));
end
