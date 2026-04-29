function [pxyz,pind,pvalue]=patch_plane4(Nodes,NodeNum,Values)

%PATCH_PLANE4  Patch information of the plane4 elements for plotting.
%
%   [pxyz,pind,pvalue] = patch_plane4(Nodes,NodeNum,Values) returns matrices
%   to plot patches of plane4 elements.
%
%   Nodes      Node definitions        [NodID x y z]
%   NodeNum    Node numbers       [NodID1 NodID2 NodID3 NodID4] (nElem * 4)
%   Values     Values assigned to nodes used for coloring    (nElem * 4)
%   pxyz       Coordinates of Nodes                          (4*nElem * 3)
%   pind       indices of Nodes                                (nElem * 4)
%   pvalue     Values arranged per Node                      (4*nElem * 1)
%
%   See also PLOTSTRESSCONTOURF, PLOTSHELLFCONTOURF.

% Miche Jansen
% 2010

[pxyz,pind,pvalue]=patch_shell4(Nodes,NodeNum,Values);

end