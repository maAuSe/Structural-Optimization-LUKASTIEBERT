function DLoads=accel_solid8(Accelxyz,Nodes,Elements,Sections,Materials,Options)

%ACCEL_SOLID8   Compute the distributed loads for solid8 elements due to an
%               acceleration.
%
%   DLoads = accel_solid8(Accelxyz,[],Elements,Sections,Materials,Options)
%   computes the distributed loads for solid8 elements due to an acceleration.
%   In order to simulate gravity, accelerate the structure in the direction
%   opposite to gravity.
%
%   Accelxyz   Acceleration            [Ax Ay Az] (1 * 3)
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Sections   Section definitions     [SecID SecProp1 SecProp2 ...]
%   Materials  Material definitions    [MatID MatProp1 MatProp2 ... ]
%   Options    Element options            {Option1 Option2 ...}
%   DLoads     Distributed loads       [EltID n1globalX n1globalY n1globalZ ...]
%
%   See also ACCEL, ACCEL_TRUSS.


% PREPROCESSING
Accelxyz=Accelxyz(:).';
nElem=size(Elements,1);

DLoads=zeros(nElem,25);
rhoInd = 4;

for iElem=1:nElem
    MatID=Elements(iElem,4);
    loc=find(Materials(:,1)==MatID);
    if isempty(loc)
        error('Material %i is not defined.',MatID)
    elseif length(loc)>1
        error('Material %i is multiply defined.',MatID)
    end
    rho=Materials(loc,rhoInd);

    DLoads(iElem,:)=[Elements(iElem,1) -rho*kron(ones(1,8),Accelxyz)];
end