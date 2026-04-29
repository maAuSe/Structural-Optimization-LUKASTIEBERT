function [S,dSdx] = size_plane6(Node,dNodedx)

% SIZE_PLANE6   Compute plane6 element size (area).
%
%   S = SIZE_PLANE6(Node) computes the size (area) of a plane6 element.
%
%   [S,dSdx] = SIZE_PLANE6(Node,dNodedx) additionally computes the
%   derivatives of the element size with respect to the design variables x.
%
%   Node        Node definitions               	[x y z] (3 * 3)
%   dNodedx     Node definitions derivatives   	(SIZE(Node) * nVar)
%   S           Element size
%   dSdx        Element size derivatives
%
%   See also ELEMSIZES, ELEMVOLUMES, SIZE_BEAM.

% Mattias Schevenels
% January 2024



if nargin>1
  nVar = size(dNodedx,3);
else
  nVar = 0;
end

S = 0;
dSdx = zeros(nVar,1);

[xi,H]=gaussqtri(7);

for iXi=1:7
    [Ni,dNi_dxi,dNi_deta]=sh_t6(xi(iXi,1),xi(iXi,2));
    J=[dNi_dxi dNi_deta].'*Node(1:6,1:2);
    detJ=det(J);
    S=S+H(iXi)*detJ;

    if nargout>1
      G = J\[dNi_dxi dNi_deta].';
      X = Node(1:6,1:2);
      for iVar = 1:nVar
        dXdx = dNodedx(1:6,1:2,iVar);
        ddetJdx = detJ*trace(G*dXdx);
        dSdx(iVar) = dSdx(iVar) + ddetJdx*H(iXi);
      end
    end
end

