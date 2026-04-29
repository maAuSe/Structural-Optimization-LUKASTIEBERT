function printdisp(Nodes,DOF,U)

%PRINTDISP   Display the displacements in the command window.
%
%    printdisp(Nodes,DOF,U)
%   displays the displacements in the command window.
%
%   Nodes      Node definitions          [NodID x y z]
%   DOF        Degrees of freedom  (nDOF * 1)
%   U          Displacements (nDOF * 1)
%
%   See also PRINTFORC, PLOTDISP.

% David Dooms
% March 2008

% PREPROCESSING
DOF=DOF(:);

nNode=size(Nodes,1);

NodeList=sort(Nodes(:,1));

NodeDOF=floor(DOF);
DOFDOF=round(rem(DOF,1)*100);

sdisp='   NODE           Ux           Uy           Uz           Rx           Ry           Rz';

for iNode=1:nNode
    ind=find(NodeDOF==NodeList(iNode));
    UNode=U(ind,1);
    indDOF=DOFDOF(ind);
    
    tableU=zeros(1,6);
    tableU(indDOF)=UNode;

    sdisp=sprintf('%s \n%7i %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e',sdisp,NodeList(iNode),tableU);
end

disp(sdisp);
