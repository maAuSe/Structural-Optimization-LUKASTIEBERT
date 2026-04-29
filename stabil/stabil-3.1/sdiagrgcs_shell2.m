function [ElemGCS,SdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema] = sdiagrgcs_shell2(stype,Forces,Node,Section,Material,DLoad,Points)

%SDIAGRGCS_SHELL2   Return matrices to plot the stresses in a SHELL2 element.
%
%   [ElemGCS,SdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                   = sdiagrgcs_shell2(ftype,Forces,Node,Section,Material,DLoad,Points)
%   [ElemGCS,SdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema]
%                   = sdiagrgcs_shell2(ftype,Forces,Node,Section,Material,DLoad)
%   returns the coordinates of the points along the SHELL2 in the global
%   coordinate system and the coordinates of the stresses with respect to the element
%   in the global coordinate system. These can be added in order to plot the
%   stresses: ElemGCS+SdiagrGCS. The coordinates of the points with extreme values
%   and the coordinates of the extreme values with respect to the element are given
%   as well and can be similarly added:  ElemExtGCS+ExtremaGCS. Extrema is the
%   list with the correspondig extreme values.
%
%   ftype      'sNphi'      Stress due to normal force in meridional direction
%              'sMphiT'     Stress at the top due to bending moment in 
%                                                              meridional direction
%              'sMphiB'     Stress at the bottom due to bending moment in 
%                                                              meridional direction
%              'sNtheta'    Stress due to normal force in circumferential direction
%              'sMthetaT'   Stress at the top due to bending moment in 
%                                                          circumferential direction
%              'sMthetaB'   Stress at the bottom due to bending moment in 
%                                                          circumferential direction
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

[ElemGCS,SdiagrGCS,ElemExtGCS,ExtremaGCS,Extrema] = fdiagrgcs_shell2(stype,Forces,Node,Section,Material,DLoad,Points);
