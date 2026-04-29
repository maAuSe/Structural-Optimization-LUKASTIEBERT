function DLoads=accel_shell2(Accelxyz,Nodes,Elements,Sections,Materials,Options)

%ACCEL_SHELL2   Compute the distributed loads for a SHELL2 element due to an acceleration.
%
%   DLoads=accel_shell2(Accelxyz,Nodes,Elements,Sections,Materials,Options)
%   computes the distributed loads for a SHELL2 element due to an acceleration.
%   In order to simulate gravity, accelerate the structure in the direction
%   opposite to gravity.
%
%   Accelxyz   Acceleration            [Ax Ay Az] (1 * 3)
%   Nodes      Node definitions        [NodeID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Sections   Section definitions     [SecID SecProp1 SecProp2 ...]
%   Materials  Material definitions    [MatID MatProp1 MatProp2 ... ]
%   Options    Element options            {Option1 Option2 ...}
%   DLoads     Distributed loads       [EltID n1globalX n1globalY n1globalZ ...]

% Difference with the ANSYS command ACEL applied to SHELL61 elements: in ANSYS
% the eccentricity of the gravity load caused by the curvature of the elements
% (more weight at the outer side, less weight at the inner side) is taken into
% account, in this function it isn't (for consistency with the element stiffness
% matrix, where a similar approximation is made).

% Mattias Schevenels
% April 2020

% PREPROCESSING
Accelxyz=Accelxyz(:).';

if Accelxyz(3)~=0, error('SHELL2 elements cannot carry loads due to an acceleration in the z-direction.'); end

nElem=size(Elements,1);

% DLoads = zeros(nElem,13); % use this line to account for eccentricity, as in ANSYS (see also kelcs_shell2)
DLoads = zeros(nElem,7);

for iElem=1:nElem
    Node1ID = Elements(iElem,5);
    loc = find(Nodes(:,1)==Node1ID);
    if isempty(loc)
        error('Node %i is not defined.',Node1ID)
    elseif length(loc)>1
        error('Node %i is multiply defined.',Node1ID)
    end
    Node2ID = Elements(iElem,6);
    loc = find(Nodes(:,1)==Node2ID);
    if isempty(loc)
        error('Node %i is not defined.',Node2ID)
    elseif length(loc)>1
        error('Node %i is multiply defined.',Node2ID)
    end
    Node = Nodes([Node1ID,Node2ID],2:4);

    SecID=Elements(iElem,3);
    loc=find(Sections(:,1)==SecID);
    if isempty(loc)
        error('Section %i is not defined.',SecID)
    elseif length(loc)>1
        error('Section %i is multiply defined.',SecID)
    end
    h=Sections(loc,2);

    MatID=Elements(iElem,4);
    loc=find(Materials(:,1)==MatID);
    if isempty(loc)
        error('Material %i is not defined.',MatID)
    elseif length(loc)>1
        error('Material %i is multiply defined.',MatID)
    end
    rho=Materials(loc,4);

    r1 = Node(1,1);
    r2 = Node(2,1);
    phi = atan2(Node(2,2)-Node(1,2),Node(2,1)-Node(1,1));

    px = -rho*Accelxyz(1);
    py = -rho*Accelxyz(2);
    Fx1 = 2*pi*r1*h*px;
    Fy1 = 2*pi*r1*h*py;
    Mz1 = 2*pi*h^3*sin(phi)/12*(px*cos(phi)+py*sin(phi));
    Fx2 = 2*pi*r2*h*px;
    Fy2 = 2*pi*r2*h*py;
    Mz2 = 2*pi*h^3*sin(phi)/12*(px*cos(phi)+py*sin(phi));


    % DLoads(iElem,[1,2,3,7,8,9,13]) = [Elements(iElem,1) Fx1 Fy1 Mz1 Fx2 Fy2 Mz2]; % use this line to account for eccentricity, as in ANSYS (see also kelcs_shell2)
    DLoads(iElem,[1,2,3,5,6]) = [Elements(iElem,1) Fx1 Fy1 Fx2 Fy2];


end
