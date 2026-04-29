function [SeGCS,SeLCS,vLCS] = se_plane6(Node,Section,Material,UeGCS,Options,gcs)
%SE_PLANE6   Compute the element stresses for a plane6 element.
%
%   [SeGCS,SeLCS,vLCS] = se_plane6(Node,Section,Material,UeGCS,Options,GCS)
%   [SeGCS,SeLCS]      = se_plane6(Node,Section,Material,UeGCS,Options,GCS)
%    SeGCS             = se_plane6(Node,Section,Material,UeGCS,Options,GCS)
%   computes the element stresses in the global and the
%   local coordinate system for the plane6 element.
%
%   Node       Node definitions           [x y z] (8 * 3)
%              Nodes should have the following order:
%              3
%              | \
%              6  5
%              |    \
%              1--4--2
%
%   Section    Section definition         [h]  (only used in plane stress)
%   Material   Material definition        [E nu rho]
%   UeGCS      Displacements (8 * nTimeSteps)
%   Options    Element options            {Option1 Option2 ...}
%   GCS        Global coordinate system in which stresses are returned
%              'cart'|'cyl'
%   SeGCS      Element stresses in GCS in corner nodes IJKL
%              48 = 6 stress comp. * 8 nodes (24 * nTimeSteps)
%                                        [sxx syy szz sxy syz sxz]
%   SeLCS      Element stresses in LCS in corner nodes IJKL
%              48 = 6 stress comp. * 8 nodes (24 * nTimeSteps)
%                                        [sxx syy szz sxy syz sxz]
%   vLCS       Unit vectors of LCS (1 * 9)
%
%   See also ELEMSTRESS, SE_PLANE8.

% Options
if nargin<5, Options=[]; end
if ~isfield(Options,'problem'), Options.problem='2dstress'; end

SeLCS = selcs_plane6(Node,Section,Material,UeGCS,Options);
SeGCS = SeLCS;

if ~isempty(gcs)
  switch lower(gcs)
  case  'cart'
    SeGCS=SeLCS;
  case  'cyl'
    a1 = [Node(1:4,1) Node(1:4,2) zeros(4,1)];
    a2 = [-Node(1:4,2) Node(1:4,1) zeros(4,1)];
  
     for iNode=1:4
     A= [a1(iNode,:)/norm(a1(iNode,:));a2(iNode,:)/norm(a2(iNode,:));0 0 1];
     theta = vtrans_solid(A);
        for z = 1:3
        SeGCS((1:6)+6*(iNode-1)+24*(z-1),:)= theta*SeGCS((1:6)+6*(iNode-1)+24*(z-1),:);
        end
     end
  case 'sph'
    error('Sperhical coordinate system not defined for 2D problem.');        
  end
end

vLCS=[];
end
