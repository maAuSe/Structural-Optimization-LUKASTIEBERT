function t = trans_solid8(Node,Options,varargin)

%TRANS_SOLID8   Transform coordinate system of a solid8 element.
%
%   t   = trans_solid8(Node)
%   t   = trans_solid8(Node,Options)
%   computes the transformation matrix between the local and the global
%   coordinate system for stress computations.
%
%   Node       Node definitions           [x y z]   (8 * 3)
%   t          Transformation matrix                (3 * 3)
%
%   See also SE_SOLID8, TRANS_TRUSS.

% Miche Jansen
% 2013

t = zeros(3,3);
t(1,:) = Node(2,:)-Node(1,:);
t(2,:) = Node(4,:)-Node(1,:);
t(3,:) = cross(t(1,:),t(2,:));
t(2,:) = cross(t(3,:),t(1,:));
t = diag(1./sqrt(sum(t.^2,2)))*t;
end

function z = cross(x,y)
z = x([2 1 1]).*y([3 3 2])-y([2 1 1]).*x([3 3 2]);
z(2) = -z(2);
end