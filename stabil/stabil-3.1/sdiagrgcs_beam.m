function [ElemGCS,SdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema] = sdiagrgcs_beam(stype,Forces,Node,Section,Material,DLoad,Points)

%SDIAGRGCS_BEAM   Return matrices to plot the stresses in a beam element.
%
%   [ElemGCS,SdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                               = sdiagrgcs_beam(stype,Forces,Node,Section,[],DLoad,Points)
%   [ElemGCS,SdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                               = sdiagrgcs_beam(Forces,Node,Section,[],DLoad)
%   returns the coordinates of the points along the beam in the global
%   coordinate system and the coordinates of the stresses with respect to the beam
%   in the global coordinate system. These can be added in order to plot the
%   stresses: ElemGCS+SdiagrGCS. The coordinates of the points with extreme values
%   and the coordinates of the extreme values with respect to the beam are given
%   as well and can be similarly added:  ElemExtGCS+ExtremaGCS. Extrema is the
%   list with the corresponding extreme values.
%
%   stype      'snorm'      Normal stress due to normal force
%              'smomyt'     Normal stress due to bending moment
%                           around the local y-direction at the top
%              'smomyb'     Normal stress due to bending moment
%                           around the local y-direction at the bottom
%              'smomzt'     Normal stress due to bending moment
%                           around the local z-direction at the top
%              'smomzb'     Normal stress due to bending moment
%                           around the local z-direction at the bottom
%              'smax'       Maximal normal stress (normal force and bending moment)
%              'smin'       Minimal normal stress (normal force and bending moment)
%   Forces     Element forces in LCS (beam convention) [N; Vy; Vz; T; My; Mz]
%                                                                        (12 * 1)
%   Node       Node definitions        [x y z] (3 * 3)
%   Section    Section definition      [A ky kz Ixx Iyy Izz yt yb zt zb]
%   DLoad      Distributed loads       [n1globalX; n1globalY; n1globalZ; ...]
%                                                                         (6 * 1)
%   Points     Points in the local coordinate system (1 * nPoints)
%   ElemGCS    Coordinates of the points along the beam in GCS (nPoints * 3)
%   SdiagrGCS  Coordinates of the stress with respect to the beam in GCS
%                                                                   (nValues * 3)
%   ElemExtGCS Coordinates of the points with extreme values in GCS (nValues * 3)
%   ExtremaGCS Coordinates of the extreme values with respect to the beam in GCS
%                                                                   (nValues * 3)
%   Extrema    Extreme values (nValues * 1)
%
%   See also PLOTSTRESS, SDIAGRLCS_BEAM, SDIAGRGCS_TRUSS.

% David Dooms
% November 2008

% PREPROCESSING
DLoad=DLoad(:);

if nargin<6
    nPoints=21;
    Points=linspace(0,1,nPoints);
else
    Points=Points(:).';
    nPoints=size(Points,2);
end

% transform distributed loads from global to local coordinate system
t=trans_beam(Node);
T=blkdiag(t,t);
DLoadLCS=T*DLoad;

% compute element length
L=norm(Node(2,:)-Node(1,:));

[SdiagrLCS,loc,Extrema] = sdiagrlcs_beam(stype,Forces,DLoadLCS,L,Section,Points);

% global coordinate system

temp=polyval([-1  1],Points);
ElemGCS=[temp.' (1-temp).']*Node(1:2,:);

temp2=polyval([-1  1],loc);
ElemExtGCS=[temp2 (1-temp2)]*Node(1:2,:);

switch lower(stype)
   case {'snorm','smomzt','smomzb','smax','smin'}       % plot in local y-plane
      SdiagrLCS=[zeros(nPoints,1) SdiagrLCS.' zeros(nPoints,1)];
      ExtremaLCS=[zeros(size(Extrema)) Extrema zeros(size(Extrema))];
   case {'smomyt','smomyb'}                               % plot in local z-plane
      SdiagrLCS=[zeros(nPoints,2) SdiagrLCS.'];
      ExtremaLCS=[zeros(size(Extrema)) zeros(size(Extrema)) Extrema];
end

SdiagrGCS=SdiagrLCS*t;
ExtremaGCS=ExtremaLCS*t;
