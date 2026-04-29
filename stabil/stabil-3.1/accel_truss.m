function DLoads=accel_truss(Accelxyz,Nodes,Elements,Sections,Materials,Options)

%ACCEL_TRUSS   Compute the distributed loads for a truss due to an acceleration.
%
%   DLoads=accel_truss(Accelxyz,[],Elements,Sections,Materials,Options)
%   computes the distributed loads for a truss due to an acceleration.
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
%   See also ACCEL, ACCEL_BEAM.

% David Dooms
% October 2008

% PREPROCESSING
Accelxyz=Accelxyz(:).';

nElem=size(Elements,1);

DLoads=zeros(nElem,7);

for iElem=1:nElem

    SecID=Elements(iElem,3);
    loc=find(Sections(:,1)==SecID);
    if isempty(loc)
        error('Section %i is not defined.',SecID)
    elseif length(loc)>1
        error('Section %i is multiply defined.',SecID)
    end
    A=Sections(loc,2);

    MatID=Elements(iElem,4);
    loc=find(Materials(:,1)==MatID);
    if isempty(loc)
        error('Material %i is not defined.',MatID)
    elseif length(loc)>1
        error('Material %i is multiply defined.',MatID)
    end
    rho=Materials(loc,4);

    DLoads(iElem,:)=[Elements(iElem,1) -Accelxyz*rho*A -Accelxyz*rho*A];

end
