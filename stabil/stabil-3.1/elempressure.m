function F=elempressure(Pressures,Nodes,Elements,Types,DOF)

%ELEMPRESSURE   Equivalent nodal forces for a pressure load on a shell element.
%
%   F = elempressure(Pressures,Nodes,Elements,Types,DOF)
%   computes the equivalent nodal forces of a distributed pressure 
%   (in the local coordinate system xyz, with z perpendicular to the surface).
%
%   Pressures  Pressure on surface or edge   [EltID n1localZ n2localZ n3localZ ...]
%   Nodes      Node definitions              [NodID x y z]
%   Elements   Element definitions           [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions      {TypID EltName Option1 ... }
%   DOF        Degrees of freedom            (nDOF * 1)
%   F          Load vector                   (nDOF * 1)
%
%   See also PRESSURE_SHELL8, PRESSURE_SHELL4, NODALVALUES.

% Miche Jansen
% 2009

nPressure=size(Pressures,1);
nterm=0;


for iPr=1:nPressure
    EltID=Pressures(iPr,1);
    loc=find(Elements(:,1)==EltID);
    if isempty(loc)
        error('Element %i is not defined.',EltID)
    elseif length(loc)>1
        error('Element %i is multiply defined.',EltID)
    end

    TypID=Elements(loc,2);
    NodeNum=Elements(loc,5:end);
    
    Node=zeros(length(NodeNum),3);
    for iNode=1:length(NodeNum)
        loc=find(Nodes(:,1)==NodeNum(1,iNode));
        if isempty(loc)
            Node(iNode,:)=[NaN NaN NaN];
        elseif length(loc)>1
            error('Node %i is multiply defined.',NodeNum(1,iNode))
        else
            Node(iNode,:)=Nodes(loc,2:end);
        end
    end
    
    loc=find(cell2mat(Types(:,1))==TypID);
    if isempty(loc)
        error('Element type %i is not defined.',TypID)
    elseif length(loc)>1
        error('Element type %i is multiply defined.',TypID)
    end
    
    Type=Types{loc,2};
    
    Pressure=permute(Pressures(iPr,2:end,:),[2 3 1]);
    
    Fe=eval(['pressure_' Type '(Pressure,Node)']);
    
    PLoad(nterm+1:nterm+length(Fe),:)=Fe;
    
    seldof(nterm+1:nterm+length(Fe),1)=eval(['dof_' Type '(NodeNum)']);
    
    nterm=nterm+length(Fe);
    
end

F=nodalvalues(DOF,seldof,PLoad);
