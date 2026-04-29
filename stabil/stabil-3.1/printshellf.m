function printshellf(Elements,F)

%PRINTSHELLF Display forces/moments in command window (shell elements).
%
%   printshellf(Elements,F)
%   displays shell forces in command window.
%
%   Elements   Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   F          Shellf matrix          [Nx Ny Nxy Mx My Mxy Vx Vy] (nElem * 32)
%
%   See also PRINTSTRESS.

nElem=size(Elements,1);

[ElemList,ind]=sort(Elements(:,1));
Shellf=F(ind,:);


sshelf='ELEM  Node            Nx           Ny          Nxy           Mx           My          Mxy           Vx           Vy';

for iElem=1:nElem
    sshelf=sprintf('%s \n %5d  %3d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e',sshelf,Elements(iElem,1),Elements(iElem,5),Shellf(iElem,(1:8)));
    sshelf=sprintf('%s \n      %5d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e',sshelf,Elements(iElem,6),Shellf(iElem,(9:16)));
    sshelf=sprintf('%s \n      %5d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e',sshelf,Elements(iElem,7),Shellf(iElem,(17:24)));
    if size(Elements,2)>=8
    sshelf=sprintf('%s \n      %5d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e',sshelf,Elements(iElem,8),Shellf(iElem,(25:32)));
    end
end

disp(sshelf);