function [ElemGCS,FdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema] = fdiagrgcs_truss(ftype,Forces,Node,Section,Material,DLoad,Points)

%FDIAGRGCS_TRUSS   Return matrices to plot the forces in a truss element.
%
%   [ElemGCS,FdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                              = fdiagrgcs_truss(ftype,Forces,Node,[],[],[],Points)
%   [ElemGCS,FdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                              = fdiagrgcs_truss(ftype,Forces,Node)
%   returns the coordinates of the points along the truss in the global
%   coordinate system and the coordinates of the forces with respect to the truss
%   in the global coordinate system. These can be added in order to plot the
%   forces: ElemGCS+FdiagrGCS. The coordinates of the points with extreme values
%   and the coordinates of the extreme values with respect to the truss are given
%   as well and can be similarly added:  ElemExtGCS+ExtremaGCS. Extrema is the
%   list with the correspondig extreme values.
%
%   ftype      'norm'       Normal force (in the local x-direction)
%   Forces     Element forces in LCS [N; 0; 0; 0; 0; 0](12 * 1)
%   Node       Node definitions        [x y z] (3 * 3)
%   Points     Points in the local coordinate system (1 * nPoints)
%   ElemGCS    Coordinates of the points along the truss in GCS (nPoints * 3)
%   FdiagrGCS  Coordinates of the force with respect to the truss in GCS
%                                                                   (nValues * 3)
%   ElemExtGCS Coordinates of the points with extreme values in GCS (nValues * 3)
%   ExtremaGCS Coordinates of the extreme values with respect to the truss in GCS
%                                                                   (nValues * 3)
%   Extrema    Extreme values (nValues * 1)
%
%   See also PLOTFORC, FDIAGRGCS_BEAM.

% David Dooms
% October 2008

if nargin<7
    nPoints=21;
    Points=linspace(0,1,nPoints);
else
    Points=Points(:).';
    nPoints=size(Points,2);
end

switch lower(ftype)
   case 'norm'
      % global coordinate system
      temp=polyval([-1  1],Points);
      ElemGCS=[temp.' (1-temp).']*Node(1:2,:);
      ElemExtGCS=Node(1:2,:);

      % transform displacements from global to local coordinate system
      t=trans_truss(Node);
      T=blkdiag(t,t);

      FdiagrLCS=[zeros(nPoints,1) Forces(1)*ones(nPoints,1) zeros(nPoints,1)];
      FdiagrGCS=FdiagrLCS*t;

      Extrema=[Forces(1); Forces(1)];
      ExtremaLCS=[zeros(size(Extrema)) Extrema zeros(size(Extrema))];
      ExtremaGCS=ExtremaLCS*t;
   case {'sheary','shearz','momx','momy','momz'}
      ElemGCS=[];
      ElemExtGCS=[];
      FdiagrGCS=[];
      Extrema=[];
      ExtremaGCS=[];
   otherwise
      error('Unknown element force.')
end
