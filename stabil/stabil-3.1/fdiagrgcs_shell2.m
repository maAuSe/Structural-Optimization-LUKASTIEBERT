function [ElemGCS,FdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema] = fdiagrgcs_shell2(ftype,Forces,Node,Section,Material,DLoad,Points)

%FDIAGRGCS_SHELL2   Return matrices to plot the forces in a SHELL2 element.
%
%   [ElemGCS,FdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                    = fdiagrgcs_shell2(ftype,Forces,Node,Section,Material,DLoad,Points)
%   [ElemGCS,FdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                    = fdiagrgcs_shell2(ftype,Forces,Node,Section,Material,DLoad)
%   returns the coordinates of the points along the SHELL2 in the global
%   coordinate system and the coordinates of the forces with respect to the element
%   in the global coordinate system. These can be added in order to plot the
%   forces: ElemGCS+FdiagrGCS. The coordinates of the points with extreme values
%   and the coordinates of the extreme values with respect to the element are given
%   as well and can be similarly added:  ElemExtGCS+ExtremaGCS. Extrema is the
%   list with the correspondig extreme values.
%
%   ftype      'Nphi'       Normal force (per unit length) in meridional direction
%              'Qphi'       Transverse force (per unit length) in meridional direction
%              'Mphi'       Bending moment (per unit length) in meridional direction
%              'Ntheta'     Normal force (per unit length) in circumferential direction
%              'Mtheta'     Bending moment (per unit length) in circumferential direction
%   Forces     Element forces in LCS (beam convention) [N; Vy; 0; 0; 0; Mz] (12 * 1)
%   Node       Node definitions        [x y z] (3 * 3)
%   DLoad      Distributed loads       [n1globalX; n1globalY; n1globalZ; ...] (6 * 1)
%   Points     Points in the local coordinate system (1 * nPoints)
%   ElemGCS    Coordinates of the points along the element in GCS (nPoints * 3)
%   FdiagrGCS  Coordinates of the force with respect to the element in GCS (nValues * 3)
%   ElemExtGCS Coordinates of the points with extreme values in GCS (nValues * 3)
%   ExtremaGCS Coordinates of the extreme values with respect to the element in GCS
%                                                                   (nValues * 3)
%   Extrema    Extreme values (nValues * 1)

% Mattias Schevenels
% April 2020

% PREPROCESSING
DLoad=DLoad(:);
if nargin<7
    nPoints=21;
    Points=linspace(0,1,nPoints);
else
    Points=Points(:).';
    nPoints=size(Points,2);
end

% TRANSFORM DISTRIBUTED LOADS FROM GLOBAL TO LOCAL COORDINATE SYSTEM
t=trans_shell2(Node);
if size(DLoad,1)==12 % distributed forces and moments
  T=blkdiag(t,t,t,t);
  DLoadLCS=T*DLoad(1:12,:);
else % only distributed forces
  T=blkdiag(t,t);
  DLoadLCS=T*DLoad(1:6,:);
end

% ELEMENT LENGTH, LOCATION, SLOPE
L = norm(Node(2,:)-Node(1,:));
r1 = Node(1,1);
phi = atan2(Node(2,2)-Node(1,2),Node(2,1)-Node(1,1));

% MATERIAL AND SECTION PROPERTIES
E = Material(1);
nu = Material(2);
h = Section(1);

[FdiagrLCS,loc,Extrema] = fdiagrlcs_shell2(ftype,Forces,DLoadLCS,r1,phi,L,h,E,nu,Points);

% GLOBAL COORDINATE SYSTEM
temp=polyval([-1  1],Points);
ElemGCS=[temp.' (1-temp).']*Node(1:2,:);

temp2=polyval([-1  1],loc);
ElemExtGCS=[temp2 (1-temp2)]*Node(1:2,:);

FdiagrLCS=[zeros(nPoints,1) FdiagrLCS.' zeros(nPoints,1)];
ExtremaLCS=[zeros(size(Extrema)) Extrema zeros(size(Extrema))];

FdiagrGCS=FdiagrLCS*t;
ExtremaGCS=ExtremaLCS*t;
