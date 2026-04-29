function printforc(Elements,Forces)

%PRINTFORC   Display the forces in the command window.
%
%    printforc(Elements,Forces)
%   displays the forces in the command window.
%
%   Elements   Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Forces     Element forces          [N Vy Vz T My Mz] (nElem * 12)
%
%   See also PRINTDISP.

% David Dooms
% March 2008

sforces=' ELEM              N           Vy           Vz            T           My           Mz';

nElem=size(Elements,1);

% (ElemList,ind)=sort(Elements(:,1));
% Forces=Forces(ind);

for iElem=1:nElem
    sforces=sprintf('%s \n%5i I %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e',sforces,Elements(iElem,1),Forces(iElem,1:6));
    sforces=sprintf('%s \n      J %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e',sforces,Forces(iElem,7:12));
end

disp(sforces);
