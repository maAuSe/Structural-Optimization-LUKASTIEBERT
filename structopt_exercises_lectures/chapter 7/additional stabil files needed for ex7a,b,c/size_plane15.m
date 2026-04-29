function [S,dSdx] = size_plane15(Node,dNodedx)

% SIZE_PLANE15   Compute plane15 element size (area).
%
%   S = SIZE_PLANE15(Node) computes the size (area) of a plane15 element.
%
%   [S,dSdx] = SIZE_PLANE15(Node,dNodedx) additionally computes the
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



if nargin>1
  nVar = size(dNodedx,3);
else
  nVar = 0;
end

S = 0;
dSdx = zeros(nVar,1);

[xi,H]=gaussqtri(13);

for iXi=1:13
    [Ni,dNi_dxi,dNi_deta]=sh_t15(xi(iXi,1),xi(iXi,2));
    J=[dNi_dxi dNi_deta].'*Node(1:15,1:2);
    detJ=det(J);
    S=S+H(iXi)*detJ;

    if nargout>1
      G = J\[dNi_dxi dNi_deta].';
      X = Node(1:15,1:2);
      for iVar = 1:nVar
        dXdx = dNodedx(1:15,1:2,iVar);
        ddetJdx = detJ*trace(G*dXdx);
        dSdx(iVar) = dSdx(iVar) + ddetJdx*H(iXi);
      end
    end
end

