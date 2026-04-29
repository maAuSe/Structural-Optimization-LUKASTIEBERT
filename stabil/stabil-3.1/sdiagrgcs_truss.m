function [ElemGCS,SdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema] = sdiagrgcs_truss(stype,Forces,Node,Section,Material,DLoad,Points)

%SDIAGRGCS_TRUSS   Return matrices to plot the stresses in a truss element.
%
%   [ElemGCS,SdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                              = sdiagrgcs_truss(stype,Forces,Node,Section,[],[],Points)
%   [ElemGCS,SdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                              = sdiagrgcs_truss(stype,Forces,Node,Section)
%   returns the coordinates of the points along the truss in the global
%   coordinate system and the coordinates of the stresses with respect to the truss
%   in the global coordinate system. These can be added in order to plot the
%   stresses: ElemGCS+SdiagrGCS. The coordinates of the points with extreme values
%   and the coordinates of the extreme values with respect to the truss are given
%   as well and can be similarly added:  ElemExtGCS+ExtremaGCS. Extrema is the
%   list with the corresponding extreme values.
%
%   stype      'snorm'      Normal stress due to normal force
%   Forces     Element forces in LCS [N; 0; 0; 0; 0; 0](12 * 1)
%   Node       Node definitions        [x y z] (3 * 3)
%   Section    Section definition      [A ky kz Ixx Iyy Izz yt yb zt zb]
%   Points     Points in the local coordinate system (1 * nPoints)
%   ElemGCS    Coordinates of the points along the truss in GCS (nPoints * 3)
%   SdiagrGCS  Coordinates of the force with respect to the truss in GCS
%                                                                   (nValues * 3)
%   ElemExtGCS Coordinates of the points with extreme values in GCS (nValues * 3)
%   ExtremaGCS Coordinates of the extreme values with respect to the truss in GCS
%                                                                   (nValues * 3)
%   Extrema    Extreme values (nValues * 1)
%
%   See also PLOTSTRESS, SDIAGRGCS_BEAM.

% David Dooms
% November 2008

if nargin<7
    nPoints=21;
    Points=linspace(0,1,nPoints);
else
    Points=Points(:).';
    nPoints=size(Points,2);
end

switch lower(stype)
   case {'snorm','smax','smin'}
      % global coordinate system
      temp=polyval([-1  1],Points);
      ElemGCS=[temp.' (1-temp).']*Node(1:2,:);
      ElemExtGCS=Node(1:2,:);

      % transform displacements from global to local coordinate system
      t=trans_truss(Node);
      T=blkdiag(t,t);

      SdiagrLCS=[zeros(nPoints,1) Forces(1)/Section(1)*ones(nPoints,1) zeros(nPoints,1)];
      SdiagrGCS=SdiagrLCS*t;

      Extrema=[Forces(1)/Section(1); Forces(1)/Section(1)];
      ExtremaLCS=[zeros(size(Extrema)) Extrema zeros(size(Extrema))];
      ExtremaGCS=ExtremaLCS*t;
   case {'smomyt','smomyb','smomzt','smomzb'}
      ElemGCS=[];
      ElemExtGCS=[];
      SdiagrGCS=[];
      Extrema=[];
      ExtremaGCS=[];
   otherwise
      error('Unknown stress.')
end
