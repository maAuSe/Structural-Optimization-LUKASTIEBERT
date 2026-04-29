function h=plotnodes(Nodes,varargin)

%PLOTNODES   Plot the nodes.
%
%    plotnodes(Nodes)
%   plots the nodes.
%
%   Nodes      Node definitions        [NodID x y z]
%
%   plotnodes(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'Numbering'    Plots the node numbers.  Default: 'on'.
%   'GCS'          Plots the global coordinate system.  Default: 'on'.
%   'Handle'       Plots in the axis with this handle.  Default: current axis.
%   Additional parameters are redirected to the PLOT3 function which plots 
%   the nodes.
%
%   h = PLOTNODES(...) returns a struct h with handles to all the objects
%   in the plot.
%
%   See also PLOTELEM, PLOTDISP.

% David Dooms
% March 2008

% PARAMETERS
if nargin<2
    paramlist={'k.'};
elseif rem(length(varargin),2)==0
    paramlist={'k.' varargin{:}};
else
    paramlist=varargin;
end
[Numbering,paramlist]=cutparam('Numbering','on',paramlist);
[GCS,paramlist]=cutparam('GCS','on',paramlist);
[haxis,paramlist]=cutparam('Handle',gca,paramlist);

% PREPARE FIGURE
nextplot=get(haxis,'NextPlot');

% PLOT NODES
h.nodes=plot3(Nodes(:,2),Nodes(:,3),Nodes(:,4),paramlist{:},'LineStyle','none','Clipping','off');

% PLOT NODE NUMBERS
if strcmpi(Numbering,'on')
    h.nodenumbers=text(Nodes(:,2),Nodes(:,3),Nodes(:,4),num2str(Nodes(:,1)),'HorizontalAlignment','right','VerticalAlignment','bottom');
    set(h.nodenumbers,'Color',get(h.nodes,'Color'));
end

% VIEWPOINT AND AXES FIGURE
if strcmpi(nextplot,'replace')
    if all(Nodes(:,4)==0)
        view(0,90);
    elseif all(Nodes(:,3)==0)
        view(0,0);
    elseif all(Nodes(:,2)==0)
        view(90,0);
    else
        view(37.5,30);
    end
    defaultposition=[0.13 0.11 0.775 0.815];
    if all(get(haxis,'Position')==defaultposition)
        set(haxis,'Position',[0.05 0.05 0.9 0.87])
    end
    axis equal
    set(gcf,'resizefcn',@figResize)
    axis off
    set(gcf,'PaperPositionMode','auto')
end

% PLOT GLOBAL COORDINATE SYSTEM
if strcmpi(GCS,'on')
    lref=reflength(Nodes);
    set(haxis,'NextPlot','add');
    h=plotgcs(lref,h);
end

% RESET NEXTPLOT STATE
set(haxis,'NextPlot',nextplot);

% RETURN OUTPUT ARGUMENTS ONLY IF REQUESTED
if nargout<1, clear('h'); end