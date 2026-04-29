function DLoads=accel_shell6(Accelxyz,Nodes,Elements,Sections,Materials,Options)

%ACCEL_SHELL6   Compute the distributed loads for shell6 elements due to an acceleration.
%
%   DLoads = accel_shell6(Accelxyz,[],Elements,Sections,Materials,Options)
%   computes the distributed loads for shell8 elements due to an acceleration.
%   In order to simulate gravity, accelerate the structure in the direction
%   opposite to gravity.
%
%   Accelxyz   Acceleration            [Ax Ay Az] (1 * 3)
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Sections   Section definitions     [SecID SecProp1 SecProp2 ...]
%   Materials  Material definitions    [MatID MatProp1 MatProp2 ... ]
%   Options    Element options struct. Fields:
%              -MatType: 'isotropic' (default) or 'orthotropic'
%   DLoads     Distributed loads       [EltID n1globalX n1globalY n1globalZ ...]
%
%   See also ACCEL, ACCEL_TRUSS.


% PREPROCESSING
if isfield(Options,'MatType')
    if strcmpi(Options.MatType,'orthotropic')
    rhoInd = 9;
    end
else
rhoInd = 4;
end

Accelxyz=Accelxyz(:).';

nElem=size(Elements,1);

DLoads=zeros(nElem,10);

for iElem=1:nElem

    SecID=Elements(iElem,3);
    loc=find(Sections(:,1)==SecID);
    if isempty(loc)
        error('Section %i is not defined.',SecID)
    elseif length(loc)>1
        error('Section %i is multiply defined.',SecID)
    end
    h=Sections(loc,2:end);

    MatID=Elements(iElem,4);
    loc=find(Materials(:,1)==MatID);
    if isempty(loc)
        error('Material %i is not defined.',MatID)
    elseif length(loc)>1
        error('Material %i is multiply defined.',MatID)
    end
    rho=Materials(loc,rhoInd);
    if length(h)<3
        DLoads(iElem,:)=[Elements(iElem,1) -rho*kron(h(1,[1 1 1]),Accelxyz)];
    elseif any(isnan(h(2:3)))
        DLoads(iElem,:)=[Elements(iElem,1) -rho*kron(h(1,[1 1 1]),Accelxyz)];
    else
    DLoads(iElem,:)=[Elements(iElem,1) -rho*kron(h(1,1:3),Accelxyz)];
    end
end