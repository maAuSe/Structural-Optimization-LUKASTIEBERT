function [ElemGCS,FdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema] = fdiagrgcs_beam(ftype,Forces,Node,Section,Material,DLoad,Points)

%FDIAGRGCS_BEAM   Return matrices to plot the forces in a beam element.
%
%   [ElemGCS,FdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                               = fdiagrgcs_beam(ftype,Forces,Node,[],[],DLoad,Points)
%   [ElemGCS,FdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                               = fdiagrgcs_beam(ftype,Forces,Node,[],[],DLoad)
%   returns the coordinates of the points along the beam in the global
%   coordinate system and the coordinates of the forces with respect to the beam
%   in the global coordinate system. These can be added in order to plot the
%   forces: ElemGCS+FdiagrGCS. The coordinates of the points with extreme values
%   and the coordinates of the extreme values with respect to the beam are given
%   as well and can be similarly added:  ElemExtGCS+ExtremaGCS. Extrema is the
%   list with the correspondig extreme values.
%
%   ftype      'norm'       Normal force (in the local x-direction)
%              'sheary'     Shear force in the local y-direction
%              'shearz'     Shear force in the local z-direction
%              'momx'       Torsional moment (around the local x-direction)
%              'momy'       Bending moment around the local y-direction
%              'momz'       Bending moment around the local z-direction
%   Forces     Element forces in LCS (beam convention) [N; Vy; Vz; T; My; Mz]
%                                                                        (12 * 1)
%   Node       Node definitions        [x y z] (3 * 3)
%   DLoad      Distributed loads       [n1globalX; n1globalY; n1globalZ; ...]
%                                                                         (6 * 1)
%   Points     Points in the local coordinate system (1 * nPoints)
%   ElemGCS    Coordinates of the points along the beam in GCS (nPoints * 3)
%   FdiagrGCS  Coordinates of the force with respect to the beam in GCS
%                                                                   (nValues * 3)
%   ElemExtGCS Coordinates of the points with extreme values in GCS (nValues * 3)
%   ExtremaGCS Coordinates of the extreme values with respect to the beam in GCS
%                                                                   (nValues * 3)
%   Extrema    Extreme values (nValues * 1)
%
%   See also PLOTFORC, FDIAGRLCS_BEAM, FDIAGRGCS_TRUSS.

% David Dooms
% October 2008

% PREPROCESSING
if nargin<7
    nPoints=21;
    Points=linspace(0,1,nPoints);
else
    Points=Points(:).';
    nPoints=size(Points,2);
end

% transform distributed loads from global to local coordinate system
t=trans_beam(Node);
T=blkdiag(t,t);
DLoadLCS=dloadgcs2lcs(T,DLoad);

% compute element length
L=norm(Node(2,:)-Node(1,:));

[FdiagrLCS,loc,Extrema] = fdiagrlcs_beam(ftype,Forces,DLoadLCS,L,Points);

% global coordinate system

temp=polyval([-1  1],Points);
ElemGCS=[temp.' (1-temp).']*Node(1:2,:);

temp2=polyval([-1  1],loc);
ElemExtGCS=[temp2 (1-temp2)]*Node(1:2,:);

switch lower(ftype)
   case {'norm','momx','sheary','momz'}                 % plot in local y-plane
      FdiagrLCS=[zeros(nPoints,1) FdiagrLCS.' zeros(nPoints,1)];
      ExtremaLCS=[zeros(size(Extrema)) Extrema zeros(size(Extrema))];
   case {'shearz','momy'}                               % plot in local z-plane
      FdiagrLCS=[zeros(nPoints,2) FdiagrLCS.'];
      ExtremaLCS=[zeros(size(Extrema)) zeros(size(Extrema)) Extrema];
end

FdiagrGCS=FdiagrLCS*t;
ExtremaGCS=ExtremaLCS*t;