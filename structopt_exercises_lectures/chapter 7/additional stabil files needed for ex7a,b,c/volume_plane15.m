function [V,dVdx] = volume_plane15(Node,Section,dNodedx,dSectiondx)

% VOLUME_PLANE15   Compute the volume of a plane15 element.
%
%   V = VOLUME_PLANE15(Node,Section) computes the volume of a plane15 element.
%
%   [V,dVdx] = VOLUME_PLANE15(Node,Section,dNodedx,dSectiondx) computes the
%   volume of a plane15 element, as well as the derivatives of the
%   volume with respect to the design variables x.
%
%   Node        Node definitions          [x y z] (15 * 3)
%   Sections    Section definitions       [SecID SecProp1 SecProp2 ...]
%   dNodedx     Node definitions derivatives  	(SIZE(Node) * nVar)
%   dSectionsdx Section definitions derivatives	(SIZE(Section) * nVar)
%   V           Element volume                  (1 * 1)
%   dVdx        Element volume derivatives      (nVar * 1)
%
%   See also ELEMVOLUMES, VOLUME_BEAM, ELEMSIZES, SIZE_TRUSS.

% Wouter Dillen
% December 2017

h = Section(1);
if nargout==1
  S = size_plane15(Node);
  V = h*S;
else
  if nargin<3, dNodedx = []; end
  if nargin<4, dSectiondx = []; end
  nVar = max(size(dNodedx,3),size(dSectiondx,3));
  if nVar==0 || isempty(dNodedx), dNodedx = zeros([size(Node),nVar]); end
  if nVar==0 || isempty(dSectiondx), dSectiondx = zeros([size(Section),nVar]); end
  dhdx = reshape(dSectiondx(1,1,:),[],1);
  [S,dSdx] = size_plane15(Node,dNodedx);
  V = h*S;
  dVdx = h.*dSdx + dhdx.*S;
end



