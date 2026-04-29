function DLoads=accel(Accelxyz,Nodes,Elements,Types,Sections,Materials)

%ACCEL   Compute the distributed loads due to an acceleration.
%
%   DLoads=accel(Accelxyz,Nodes,Elements,Types,Sections,Materials)
%   computes the distributed loads due to an acceleration.
%   In order to simulate gravity, accelerate the structure in the direction
%   opposite to gravity.
%
%   Accelxyz   Acceleration            [Ax Ay Az] (1 * 3)
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   Sections   Section definitions     [SecID SecProp1 SecProp2 ...]
%   Materials  Material definitions    [MatID MatProp1 MatProp2 ... ]
%   DLoads     Distributed loads       [EltID n1globalX n1globalY n1globalZ ...]
%
%   See also ELEMLOADS, ACCEL_BEAM, ACCEL_TRUSS.

% David Dooms
% October 2008

% Note (MS, April 2020):
% The SHELL2 element needs the node locations to compute acceleration loads;
% the syntax of this function is therefore modified in stabil 3.1:
% - new syntax: DLoads=accel(Accelxyz,Nodes,Elements,Types,Sections,Materials)
% - old syntax: DLoads=accel(Accelxyz,Elements,Types,Sections,Materials)
% If the old syntax is used, shift the input arguments for backward compatibility
if nargin==5
  Materials = Sections;
  Sections = Types;
  Types = Elements;
  Elements = Nodes;
  Nodes = [];
end

nType=size(Types,1);

DLoads=zeros(size(Elements,1),7);

nterm=0;

for iType=1:nType
    loc=find(Elements(:,2)==cell2mat(Types(iType,1)));
    if isempty(loc)

    else
        ElemType=Elements(loc,:);
        Type=Types{iType,2};
        if size(Types,2)<3
        Options={};
        else
        Options=Types{iType,3};
        end
        dloads2=eval(['accel_' Type '(Accelxyz,Nodes,ElemType,Sections,Materials,Options)']);
        if size(dloads2,2) > size(DLoads,2);
        DLoads = [DLoads,nan(size(Elements,1),size(dloads2,2)-size(DLoads,2))];
        DLoads((nterm+1):(nterm+length(loc)),:) = dloads2;
        elseif size(dloads2,2) < size(DLoads,2);
        DLoads((nterm+1):(nterm+length(loc)),:) = [dloads2,nan(length(loc),size(DLoads,2)-size(dloads2,2))];
        else
        DLoads((nterm+1):(nterm+length(loc)),:) = dloads2;
        end

        nterm=nterm+length(loc);
    end
end
