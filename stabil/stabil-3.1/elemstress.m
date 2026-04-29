function [SeGCS,SeLCS,vLCS]=elemstress(Nodes,Elements,Types,Sections,Materials,DOF,U,varargin)

%ELEMSTRESS   Compute the element stresses.
%
%   [SeGCS,SeLCS,vLCS] = elemstress(Nodes,Elements,Types,Sections,Materials,DOF,U)
%   [SeGCS,SeLCS]      = elemstress(Nodes,Elements,Types,Sections,Materials,DOF,U)
%    SeGCS             = elemstress(Nodes,Elements,Types,Sections,Materials,DOF,U)
%   computes the element stresses in the global and the local coordinate system.
%
%   Nodes      Node definitions          [NodID x y z]
%   Elements   Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   Sections   Section definitions       [SecID SecProp1 SecProp2 ...]
%   Materials  Material definitions      [MatID MatProp1 MatProp2 ... ]
%   DOF        Degrees of freedom  (nDOF * 1)
%   U          Displacements (nDOF * 1)
%   elemstress(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'GCS'          Sets the gcs in which SeGCS is computed.
%                  Default: 'cart'. Type of values:
%                  'cart': cartesian coordinate system
%                  'cyl' : cylindrical coordinate system
%                  'sph' : spherical coordinate system
%
%   SeGCS      Element stresses in GCS in corner nodes IJKL and 
%              at top/mid/bot of shell (nElem * 72)   
%              72 = 6 stresscomp. * 4 nodes * 3 locations
%                                        [sxx syy szz sxy syz sxz]
%   SeLCS      Element stresses in LCS in corner nodes IJKL and 
%              at top/mid/bot of shell (nElem * 72)   
%                                        [sxx syy szz sxy syz sxz]
%   vLCS       Unit vectors of LCS       (nElem * 9)
%
%   See also SE_SHELL8, SE_SHELL4.

% Miche Jansen
% 2010

if nargin<8                               
    paramlist={};
elseif nargin>5 && ischar(varargin{1})
    paramlist=varargin;
end
[gcs,paramlist]=cutparam('GCS','cart',paramlist);




nElem=size(Elements,1);
nTimeSteps=size(U,2);

SeLCS=zeros(nElem,72,nTimeSteps);
SeGCS=zeros(nElem,72,nTimeSteps);
if nargout > 2
    vLCS = nan(nElem,9);
end

for iElem=1:nElem

    TypID=Elements(iElem,2);
    loc=find(cell2mat(Types(:,1))==TypID);
    if isempty(loc)
        error('Element type %i is not defined.',TypID)
    elseif length(loc)>1
        error('Element type %i is multiply defined.',TypID)
    end
    
    Type=Types{loc,2};
    
    if size(Types,2)<3
        Options={};
    else
        Options=Types{loc,3};
    end
    
    SecID=Elements(iElem,3);
    loc=find(Sections(:,1)==SecID);
    if isempty(loc)
        error('Section %i is not defined.',SecID)
    elseif length(loc)>1
        error('Section %i is multiply defined.',SecID)
    end
    Section=Sections(loc,2:end);
    
    MatID=Elements(iElem,4);
    loc=find(Materials(:,1)==MatID);
    if isempty(loc)
        error('Material %i is not defined.',MatID)
    elseif length(loc)>1
        error('Material %i is multiply defined.',MatID)
    end
    Material=Materials(loc,2:end);
    
    NodeNum=Elements(iElem,5:end);
    
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
    
    UeGCS=elemdisp(Type,NodeNum,DOF,U);
    
    if nargout > 2
    [temp1,temp2,temp3] = eval(['se_' Type '(Node,Section,Material,UeGCS,Options,gcs)']);
    vLCS(iElem,:) = temp3;
    else
    [temp1,temp2] = eval(['se_' Type '(Node,Section,Material,UeGCS,Options,gcs)']);
    end
    SeGCS(iElem,1:size(temp1,1),:) = permute(temp1,[3 1 2]);
    SeLCS(iElem,1:size(temp2,1),:) = permute(temp2,[3 1 2]);

end