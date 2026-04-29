function t = trans_shell8(Node,Options,varargin)

%TRANS_SHELL8   Transform coordinate system of a shell8 element.
%
%   t   = trans_shell8(Node)
%   t   = trans_shell8(Node,Options)
%   computes the transformation matrix between the local and the global
%   coordinate system for stress computations.
%
%   Node       Node definitions           [x y z] (8 * 3)
%   t          Transformation matrix  (3 * 3)
%   Options    Element options struct. Fields:
%              -LCSType: determine the reference local element
%                        coordinate system. Values: 
%                        'element' (default) or 'global'
%
%   See also SE_SHELL8, TRANS_TRUSS.

% Miche Jansen
% 2013

if nargin > 1 & isfield(Options,'LCSType') 
    LCSType = Options.LCSType;    
else
    LCSType = 'element';
end

switch lower(LCSType)
    case 'element'
v10 = Node(6,:)-Node(8,:);
v20 = Node(7,:)-Node(5,:);
v10 = v10/norm(v10);
v30 = cross(v10,v20);v30 = v30/norm(v30);
v20 = cross(v30,v10);v20 = v20/norm(v20);
    case 'global'
v10 = [1 0 0];v20 = [0 1 0];v30 = [0 0 1];               
end

t = [v10;v20;v30];
end

function z = cross(x,y)
z = x([2 1 1]).*y([3 3 2])-y([2 1 1]).*x([3 3 2]);
z(2) = -z(2);
end