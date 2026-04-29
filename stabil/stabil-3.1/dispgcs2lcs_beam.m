function UeLCS=dispgcs2lcs_beam(UeGCS,Node)

%DISPGCS2LCS_BEAM   Transform the element displacements to the LCS for a beam.
%
%   UeLCS=dispgcs2lcs_beam(UeGCS,Node)
%   transforms the element displacements from the GCS to the LCS for a beam 
%   element.
%
%   Node       Node definitions           [x y z] (3 * 3)
%   UeGCS      Displacements in the GCS (12 * 1)
%   UeLCS      Displacements in the LCS (12 * 1)
%
%   See also DISPGCS2LCS_TRUSS.

% David Dooms
% September 2008

% transform displacements from global to local coordinate system
t=trans_beam(Node);
T=blkdiag(t,t,t,t);
UeLCS=T*UeGCS;
