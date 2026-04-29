function [S,dSdx] = size_beam(Node,dNodedx)

% SIZE_BEAM     Compute beam element size (length).
%
%   S = SIZE_BEAM(Node) computes the element size (length) of a two node
%   beam element.
%
%   [S,dSdx] = SIZE_BEAM(Node,dNodedx) additionally computes the
%   derivatives of the element size with respect to the design variables x.
%
%   Node        Node definitions               	[x y z] (3 * 3)
%   dNodedx     Node definitions derivatives	(SIZE(Node) * nVar)
%   S           Element size
%   dSdx        Element size derivatives
%
%   See also ELEMSIZES, ELEMVOLUMES, SIZE_TRUSS.

% Wouter Dillen
% December 2017

if nargin<2, dNodedx = []; end

nVar = 0;
if nargout>1 && ~isempty(dNodedx)
    nVar = size(dNodedx,3);
end

if nVar==0 || isempty(dNodedx), dNodedx = zeros([size(Node),nVar]); end
    
    
S = norm(Node(2,:)-Node(1,:)); 
dSdx = (permute(dNodedx(2,:,:)-dNodedx(1,:,:),[3 2 1])*(Node(2,:)-Node(1,:)).')/S;
