function [V,dVdx] = volume_truss(Node,Section,dNodedx,dSectiondx)

% VOLUME_TRUSS     Compute the volume of a truss element.
%
%   V = VOLUME_TRUSS(Node,Section) computes the volume of a two-node truss
%   element.
%
%   [V,dVdx] = VOLUME_TRUSS(Node,Section,dNodedx,dSectiondx) computes the
%   volume of a two node truss element, as well as the derivatives of the
%   volume with respect to the design variables x.
%
%   Node        Node definitions          [x y z] (3 * 3)
%   Sections    Section definitions       [SecID SecProp1 SecProp2 ...]
%   dNodedx     Node definitions derivatives  	(SIZE(Node) * nVar)
%   dSectionsdx Section definitions derivatives	(SIZE(Section) * nVar)
%   V           Element volume                  (1 * 1)
%   dVdx        Element volume derivatives      (nVar * 1)
%
%   See also ELEMVOLUMES, VOLUME_BEAM, ELEMSIZES, SIZE_TRUSS.

% Wouter Dillen
% December 2017

if nargin<3, dNodedx = []; end
if nargin<4, dSectiondx = []; end
    
[V,dVdx] = volume_beam(Node,Section,dNodedx,dSectiondx);
    
