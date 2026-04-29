function [h,DispScal] = plotdispx(Nodes,Elements,Types,DOF,U,x,varargin);
%PLOTDISPX   Plot element quantity on deformed elements.
%
%   plotdispx(Nodes,Elements,Types,DOF,U,x);
%   [h,DispScal] = plotdispx(Nodes,Elements,Types,DOF,U,x,varargin)
%    plots element quantity on deformed elements.
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   DOF        Degrees of freedom  (nDOF * 1)
%   U          Displacements (nDOF * 1)
%   x          Elements quantity to plot   (nElem * 1)
%   plotdispx(...,ParamName,ParamValue) sets the value of the 
%   specified parameters.  The following parameters can be specified:
%   'GCS'          Plot the GCS. Default: 'on'.
%   'minmax'       Add location of min and max stress. Default: 'off'.
%   'colorbar'     Add colorbar. Default: 'off'.
%   'Undeformed'   Plots the undeformed mesh.  {'on' | 'off' (default)}
%   'ncolor'       Number of colors in colormap. Default: 10.
%   'Handle'       Plots in the axis with this handle.  Default: current axis.
%   Additional parameters are redirected to the PATCH function which plots
%   the elements.

% Miche Jansen
% 2013

% PARAMETERS
[haxis,varargin]=cutparam('Handle',gca,varargin);
[DispScal,varargin]=cutparam('DispScal','auto',varargin);
[Undeformed,varargin]=cutparam('Undeformed','off',varargin);


% PREPARE FIGURE
nextplot=get(haxis,'NextPlot');

% PLOT DEFORMED MESH
Ux = selectdof(DOF,Nodes(:,1)+0.01)*U;
Uy = selectdof(DOF,Nodes(:,1)+0.02)*U;
Uz = selectdof(DOF,Nodes(:,1)+0.03)*U;

if ischar(x)
    switch lower(x)
       case 'ux'
       edof = 0.01;
       case 'uy'   
       edof = 0.02;          
       case 'uz'
       edof = 0.03;              
    end
    edof = Elements(:,5:end)+edof;
    [lue,iue] = ismember(edof,DOF);
    x = zeros(size(Elements,1),size(Elements,2)-4);
    x(lue) = U(iue(lue)); 
end

% CHECK SOLUTION
lref=reflength(Nodes);
DispAmpl=sqrt(Ux.^2+Uy.^2+Uz.^2);
MaxDisp=max(max(DispAmpl));
% AUTOMATIC DISPLACEMENT SCALING
if strcmpi(DispScal,'auto')
    if MaxDisp > lref*1e-9
    DispScal=3*lref/MaxDisp;
    else
    DispScal=0;    
    end
end

Nodesu = Nodes;
Nodesu(:,2:4) = Nodesu(:,2:4)+DispScal*[Ux,Uy,Uz];

h.disp = plotelemx(Nodesu,Elements,Types,x,'Handle',haxis,varargin{:});

% PLOT UNDEFORMED MESH
if ~strcmpi(Undeformed,'off')
    set(haxis,'NextPlot','add');
    h=plotelemx(Nodes,Elements,Types,zeros(size(Elements,1)),'FaceColor','none','LineStyle',':','Handle',haxis,'GCS','off');
end

% RESET NEXTPLOT STATE
set(haxis,'NextPlot',nextplot);

end