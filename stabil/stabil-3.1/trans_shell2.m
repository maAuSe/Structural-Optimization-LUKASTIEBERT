function t = trans_beam(Node,varargin)

%TRANS_BEAM   Transform coordinate system for a beam element.
%
%   t = trans_beam(Node)
%   computes the transformation matrix between the local and the global
%   coordinate system for the SHELL2 element.
%
%   Node       Node definitions           [x y z] (3 * 3)
%   t          Transformation matrix  (3 * 3)
%
%   See also KE_BEAM, TRANS_TRUSS.

% Mattias Schevenels
% April 2020

phi = atan2(Node(2,2)-Node(1,2),Node(2,1)-Node(1,1));
t = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];

