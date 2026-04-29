function [X,Y,Z]=coord_kbeam(Nodes,NodeNum)

%COORD_BEAM  Coordinates of the beam elements for plotting.
%
%   coord_beam(Nodes,NodeNum)
%   returns the coordinates of the beam elements for plotting.
%
%   Nodes      Node definitions        [NodID x y z] (nNodes * 4)
%   NodeNum Node numbers          [NodID1 NodID2 NodID3] (nElem * 3)
%   X          X coordinates  (2 * nElem)    
%   Y          Y coordinates  (2 * nElem)
%   Z          Z coordinates  (2 * nElem)
%
%   See also COORD_TRUSS, PLOTELEM.

% David Dooms
% March 2008

[X,Y,Z] = coord_beam(Nodes,NodeNum);
