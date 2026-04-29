function [V,dVdx] = volume_beam(Node,Section,dNodedx,dSectiondx)

% VOLUME_BEAM     Compute the volume of a beam element.
%
%   V = VOLUME_BEAM(Node,Section) computes the volume of a two-node beam
%   element.
%
%   [V,dVdx] = VOLUME_BEAM(Node,Section,dNodedx,dSectiondx) computes the
%   volume of a two node beam element, as well as the derivatives of the
%   volume with respect to the design variables x.
%
%   Node        Node definitions          [x y z] (3 * 3)
%   Sections    Section definitions       [SecID SecProp1 SecProp2 ...]
%   dNodedx     Node definitions derivatives   	(SIZE(Node) * nVar)
%   dSectionsdx Section definitions derivatives	(SIZE(Section) * nVar)
%   V           Element volume                  (1 * 1)
%   dVdx        Element volume derivatives      (nVar * 1)
%
%   See also ELEMVOLUMES, VOLUME_TRUSS, ELEMSIZES, SIZE_BEAM.

% Wouter Dillen
% December 2017

% preprocessing
if nargin<3, dNodedx = []; end
if nargin<4, dSectiondx = []; end

nVar = 0;
if nargout>1 && (~isempty(dNodedx) || ~isempty(dSectiondx))
    nVar = max(size(dNodedx,3),size(dSectiondx,3));
end

if nVar==0 || isempty(dNodedx), dNodedx = zeros([size(Node),nVar]); end
if nVar==0 || isempty(dSectiondx), dSectiondx = zeros([size(Section),nVar]); end
    

% Length
[L,dLdx] = size_beam(Node,dNodedx);

% Section area
A = Section(1);
dAdx = permute(dSectiondx(:,1,:),[3 2 1]);

% Volume
V = A*L;
dVdx = dAdx*L + A*dLdx;

