function [DOF,elemInd] = getdof(Elements,Types)

%GETDOF   Get the vector with the degrees of freedom of the model.
%
%   DOF=getdof(Elements,Types) builds the vector with the
%   labels of the degrees of freedom for which stiffness is present in the
%   finite element model.
%
%   Elements   Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   DOF        Degrees of freedom  (nDOF * 1)
%
%   See also DOF_TRUSS, DOF_BEAM, GETDOF.


nElem=size(Elements,1);

elemInd = cell(nElem,1);
DOF = cell(nElem,1);
nterm=0;

checkunique(cell2mat(Types(:,1)),'Element type');
Elements(:,2)=linkandcheck(round(Elements(:,2)),round(cell2mat(Types(:,1))),'Element type');
% Element dof functions handles
doffunc = cellfun(@(x) ['dof_' x],Types(:,2),'UniformOutput',0);
doffunc = cellfun(@str2func,doffunc,'UniformOutput',0);

for iElem=1:nElem
    
    NodeNum=Elements(iElem,5:end);
    dofelem = doffunc{Elements(iElem,2)}(NodeNum);
    
    elemInd{iElem}=nterm+(1:length(dofelem));
    DOF{iElem}=dofelem;
    nterm=nterm+length(dofelem);
end

DOF = cell2mat(DOF);
[DOF,iDOF,jDOF]=unique(DOF);

for iElem=1:nElem
    elemInd{iElem}=jDOF(elemInd{iElem});
end
end