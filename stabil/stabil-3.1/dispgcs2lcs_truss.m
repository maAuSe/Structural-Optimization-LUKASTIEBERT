function UeLCS=dispgcs2lcs_truss(UeGCS,Node)

%DISPGCS2LCS_TRUSS  Transform the element displacements to the LCS for a truss.
%
%   UeLCS=dispgcs2lcs_truss(UeGCS,Node)
%   transforms the element displacements from the GCS to the LCS for a truss 
%   element.
%
%   Node       Node definitions           [x y z] (2 * 3)
%   UeGCS      Displacements in the GCS (6 * 1)
%   UeLCS      Displacements in the LCS (6 * 1)
%
%   See also DISPGCS2LCS_BEAM.

% David Dooms
% September 2008

% transform displacements from global to local coordinate system
t=trans_truss(Node);
T=blkdiag(t,t);
UeLCS=T*UeGCS;
