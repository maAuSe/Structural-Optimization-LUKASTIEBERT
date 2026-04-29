function h=plotelem(Nodes,Elements,Types,varargin)

%PLOTELEM   Plot the elements.
%
%    plotelem(Nodes,Elements,Types)
%   plots the elements.
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%
%   plotelem(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'Numbering'    Plots the element numbers.  Default: 'on'.
%   'GCS'          Plots the global coordinate system.  Default: 'on'.
%   'Handle'       Plots in the axis with this handle.  Default: current axis.
%   Additional parameters are redirected to the PLOT3 function which plots
%   the elements.
%
%   h = PLOTELEM(...) returns a struct h with handles to all the objects
%   in the plot.
%
%   See also COORD_TRUSS, COORD_BEAM, PLOTDISP.

% David Dooms
% March 2008

% PARAMETERS
if nargin<4
    paramlist={'k'};
elseif rem(length(varargin),2)==0
    paramlist={'k' varargin{:}};
else
    paramlist=varargin;
end
[Numbering,paramlist]=cutparam('Numbering','on',paramlist);
[GCS,paramlist]=cutparam('GCS','on',paramlist);
[haxis,paramlist]=cutparam('Handle',gca,paramlist);

nType=size(Types,1);

nterm=0;
ntermnum=0;

for iType=1:nType
    TypID=Types{iType,1};
    Type=Types{iType,2};
    ind=find(Elements(:,2)==TypID);
    NodeNumbers=Elements(ind,5:end);
    ElemNumbers(ntermnum+1:ntermnum+length(ind),:)=Elements(ind,1);
    [x,y,z]=eval(['coord_' Type '(Nodes,NodeNumbers)']); % aanpassen!!!
    cx(ntermnum+1:ntermnum+length(ind))=mean(x);
    cy(ntermnum+1:ntermnum+length(ind))=mean(y);
    cz(ntermnum+1:ntermnum+length(ind))=mean(z);
    ntermnum=ntermnum+length(ind);
    if size(x,1)>2
      x=[x;x(1,:)];
      y=[y;y(1,:)];
      z=[z;z(1,:)];
    end
    x=[x; nan(1,size(x,2))];
    y=[y; nan(1,size(y,2))];
    z=[z; nan(1,size(z,2))];
    X(nterm+1:nterm+prod(size(x)),1)=x(:);
    Y(nterm+1:nterm+prod(size(x)),1)=y(:);
    Z(nterm+1:nterm+prod(size(x)),1)=z(:);
    nterm=nterm+prod(size(x));
end

% PREPARE FIGURE
nextplot=get(haxis,'NextPlot');

% PLOT ELEMENTS
h.elements=plot3(haxis,X,Y,Z,paramlist{:},'Clipping','off');


% PLOT ELEMENT NUMBERS
if strcmpi(Numbering,'on')
    h.elementnumbers=text(cx,cy,cz,num2str(ElemNumbers),'HorizontalAlignment','center','VerticalAlignment','mid','Parent',haxis);
    set(h.elementnumbers,'Color',get(h.elements,'Color'));
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