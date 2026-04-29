function [S,dSdx] = size_plane3(Node,dNodedx)

% SIZE_PLANE3   Compute plane3 element size (area).
%
%   S = SIZE_PLANE3(Node) computes the size (area) of a plane3 element.
%
%   [S,dSdx] = SIZE_PLANE3(Node,dNodedx) additionally computes the
%   derivatives of the element size with respect to the design variables x.
%
%   Node        Node definitions               	[x y z] (3 * 3)
%   dNodedx     Node definitions derivatives   	(SIZE(Node) * nVar)
%   S           Element size
%   dSdx        Element size derivatives
%
%   See also ELEMSIZES, ELEMVOLUMES, SIZE_BEAM.

% Mattias Schevenels
% November 2023

% Element area
X=Node(:,1);
Y=Node(:,2);

XY = [1 X(1) Y(1);
      1 X(2) Y(2);
      1 X(3) Y(3)];

S = 0.5*det(XY);

if nargin>1
  nVar = size(dNodedx,3);
  nul = zeros(1,1,nVar);
  dXdx = dNodedx(:,1,:);
  dYdx = dNodedx(:,2,:);
  dXYdx = [nul dXdx(1,1,:) dYdx(1,1,:);
           nul dXdx(2,1,:) dYdx(2,1,:);
           nul dXdx(3,1,:) dYdx(3,1,:)];
  dSdx = zeros(nVar,1);
  for iVar = 1:nVar
    dSdx(iVar) = 0.5*trace(adjoint(XY)*dXYdx(:,:,iVar));
  end
end
