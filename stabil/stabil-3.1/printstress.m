function printstress(Elements,S,location)

%PRINTSTRESS Display stress in command window (shell elements).
%
%   printstress(Elements,S,location) displays stress in command window.
%
%   Elements   Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   S          Stress matrix          [sx sy sz sxy syz szx] (nElem * 72)
%
%   See also PRINTSHELLF.

sstress='';
nElem=size(Elements,1);

[ElemList,ind]=sort(Elements(:,1));
Stress=S(ind,:);
Elements = Elements(ind,:);
loc={'TOP','MIDDLE','BOTTOM'};

if nargin > 2
switch lower(location)
    case {'all'}
        ind = [1,2,3];
    case {'top'}
        ind = 1;
    case {'mid'}
        ind = 2;
    case {'bot'}
        ind = 3;
end
else ind = (1:3);
end

for ind = ind 
    sstress=sprintf('%s \n LOCATION: %s',sstress,loc{ind});
    sstress=sprintf('%s \n ELEM  Node          sxx          syy          szz          sxy          syz          szx',sstress);
for iElem=1:nElem
    sstress=sprintf('%s \n %5d  %3d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e',sstress,Elements(iElem,1),Elements(iElem,5),Stress(iElem,(1:6)+24*(ind-1)));
    sstress=sprintf('%s \n      %5d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e',sstress,Elements(iElem,6),Stress(iElem,(7:12)+24*(ind-1)));
    sstress=sprintf('%s \n      %5d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e',sstress,Elements(iElem,7),Stress(iElem,(13:18)+24*(ind-1)));
    if size(Elements,2)>=8
    sstress=sprintf('%s \n      %5d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e',sstress,Elements(iElem,8),Stress(iElem,(19:24)+24*(ind-1)));
    end
end
end
disp(sstress);

