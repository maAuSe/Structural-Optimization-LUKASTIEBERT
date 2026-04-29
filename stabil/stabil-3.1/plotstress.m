function [StressScal,h]=plotstress(stype,Nodes,Elements,Types,Sections,varargin)

%PLOTSTRESS   Plot the stresses.
%
%   plotstress(stype,Nodes,Elements,Types,Sections,Materials,Forces,DLoads)
%   plotstress(stype,Nodes,Elements,Types,Sections,Materials,Forces)
%   plots the stresses.
%
%   stype      'snorm'      Normal stress due to normal force
%   (for truss 'smomyt'     Normal stress due to bending moment
%    and beam               around the local y-direction at the top
%    elements) 'smomyb'     Normal stress due to bending moment
%                           around the local y-direction at the bottom
%              'smomzt'     Normal stress due to bending moment
%                           around the local z-direction at the top
%              'smomzb'     Normal stress due to bending moment
%                           around the local z-direction at the bottom
%              'smax'       Maximal normal stress (normal force and bending moment)
%              'smin'       Minimal normal stress (normal force and bending moment)
%   stype      'sNphi'      Stress due to normal force in meridional direction
%   (for       'sMphiT'     Stress at the top due to bending moment in 
%                                                              meridional direction
%    shell2    'sMphiB'     Stress at the bottom due to bending moment in 
%                                                              meridional direction
%    elements) 'sNtheta'    Stress due to normal force in circumferential direction
%              'sMthetaT'   Stress at the top due to bending moment in 
%                                                         circumferential direction
%              'sMthetaB'   Stress at the bottom due to bending moment in 
%                                                         circumferential direction
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   Sections   Section definitions      [A ky kz Ixx Iyy Izz yt yb zt zb]
%   Materials  Material definitions      [MatID MatProp1 MatProp2 ... ]
%   Forces     Element forces in LCS (beam convention) [N Vy Vz T My Mz]
%                                                                   (nElem * 12)
%   DLoads     Distributed loads       [EltID n1globalX n1globalY n1globalZ ...]
%
%   plotstress(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'StressScal'   Stress scaling.  Default: 'auto'.
%   'MinMax'       Display minimum and maximum value.  Default: 'off'.
%   'Values'       Stress values.  Default: 'on'.
%   'Undeformed'   Plots the undeformed mesh.  Default: 'k-'.
%   'Handle'       Plots in the axis with this handle.  Default: current axis.
%   Additional parameters are redirected to the PLOT3 function which plots
%   the stresses.
%
%   [StressScal,h] = plotstress(...) returns a struct h with handles to all the
%   objects in the plot.
%
%   See also SDIAGRGCS_BEAM, SDIAGRGCS_TRUSS.

% David Dooms
% October 2008

% Note (MS, April 2020):
% The SHELL2 element needs the Materials to plot stresses along elements;
% the syntax of this function is therefore modified in stabil 3.1:
% - new syntax: plotstress(stype,Nodes,Elements,Types,Sections,Materials,Forces,DLoads,...)
% - old syntax: plotstress(stype,Nodes,Elements,Types,Sections,Forces,DLoads,...)
% If the old syntax is used, shift the input arguments for backward compatibility
if nargin==6
  newsyntax = false; % old syntax for sure
elseif isstr(varargin{2})
  newsyntax = false; % old syntax for sure
elseif nargin>=8 && ~isstr(varargin{2}) && ~isstr(varargin{3}) %
  newsyntax = true;  % new syntax for sure
else % not sure; use dimensions of Forces|Materials to differentiate (~= foolproof...)
  newsyntax = size(varargin{1},2)~=12;
end
if newsyntax
  Materials = varargin{1};
  Forces = varargin{2};
  varargin = varargin(3:end);
else
  Forces = varargin{1};
  varargin = varargin(2:end);
end

if ~isempty(varargin) && ~isstr(varargin{1})
  DLoads = varargin{1};
  varargin = varargin(2:end);
end
paramlist = varargin;

% PARAMETERS
[StressScal,paramlist]=cutparam('StressScal','auto',paramlist);
[MinMax,paramlist]=cutparam('MinMax','off',paramlist);
[Values,paramlist]=cutparam('Values','on',paramlist);
[Undeformed,paramlist]=cutparam('Undeformed','k-',paramlist);
[haxis,paramlist]=cutparam('Handle',gca,paramlist);

ElemGCS=[];
SdiagrGCS=[];
ElemExtGCS=[];
ExtremaGCS=[];
Extrema=[];

nElem=size(Elements,1);

nPoint = max(2,min(21,round(300/nElem)));
Points = linspace(0,1,nPoint);

for iElem=1:nElem

    TypID=Elements(iElem,2);
    loc=find(cell2mat(Types(:,1))==TypID);
    if isempty(loc)
        error('Element type %i is not defined.',TypID)
    elseif length(loc)>1
        error('Element type %i is multiply defined.',TypID)
    end

    Type=Types{loc,2};

    SecID=Elements(iElem,3);
    loc=find(Sections(:,1)==SecID);
    if isempty(loc)
        error('Section %i is not defined.',SecID)
    elseif length(loc)>1
        error('Section %i is multiply defined.',SecID)
    end
    Section=Sections(loc,2:end);

    if exist('Materials','var')
      MatID=Elements(iElem,4);
      loc=find(Materials(:,1)==MatID);
      if isempty(loc)
          error('Material %i is not defined.',MatID)
      elseif length(loc)>1
          error('Material %i is multiply defined.',MatID)
      end
      Material=Materials(loc,2:end);
    else
      Material=[];
    end

    NodeNum=Elements(iElem,5:end);

    Node=zeros(length(NodeNum),3);
    for iNode=1:length(NodeNum)
        loc=find(Nodes(:,1)==NodeNum(1,iNode));
        if isempty(loc)
            Node(iNode,:)=[NaN NaN NaN];
        elseif length(loc)>1
            error('Node %i is multiply defined.',NodeNum(1,iNode))
        else
            Node(iNode,:)=Nodes(loc,2:end);
        end
    end

    Force=Forces(iElem,:).';

    EltID=Elements(iElem,1);

    nTimeSteps=1;
    if ~ exist('DLoads')
        DLoads=zeros(1,7,nTimeSteps);
    else
        if size(DLoads,2)>7
            error('No partial DLoads are allowed when plotting stresses.')
        end
    end
    loc=find(DLoads(:,1)==EltID);
    if isempty(loc)
        DLoad=zeros(6,nTimeSteps);
    elseif length(loc)>1
        error('Element %i has multiple distributed loads.',EltID)
    else
        DLoad=permute(DLoads(loc,2:end,:),[2 3 1]);
    end

    [ElemGCS2,SdiagrGCS2,ElemExtGCS2,ExtremaGCS2,Extrema2]=eval(['sdiagrgcs_' Type '(stype,Force,Node,Section,Material,DLoad,Points)']);

    ElemGCS=[ElemGCS; ElemGCS2; nan(1,3)];
    SdiagrGCS=[SdiagrGCS; SdiagrGCS2; nan(1,3)];
    ElemExtGCS=[ElemExtGCS; ElemExtGCS2];
    ExtremaGCS=[ExtremaGCS; ExtremaGCS2];
    Extrema=[Extrema; Extrema2];

end

% AUTOMATIC STRESS SCALING
if strcmpi(StressScal,'auto')
    lref=reflength(Nodes);
    StressScal=3*lref/max(abs(Extrema));
end

SdiagrGCS=StressScal*SdiagrGCS+ElemGCS;
ExtremaGCS=StressScal*ExtremaGCS+ElemExtGCS;

% PREPARE FIGURE
nextplot=get(haxis,'NextPlot');

% PLOT STRESS
handlestresses=plot3(haxis,SdiagrGCS(:,1),SdiagrGCS(:,2),SdiagrGCS(:,3),paramlist{:},'Clipping','off');
set(haxis,'NextPlot','add');

% PLOT UNDEFORMED MESH
if ~ strcmpi(Undeformed,'off')
    h=plotelem(Nodes,Elements,Types,Undeformed,'Numbering','off','Handle',haxis,'GCS','off');
end
h.stresses=handlestresses;

% PLOT HATCH
h.hatching=plot3([ElemGCS(:,1).'; SdiagrGCS(:,1).'],[ElemGCS(:,2).'; SdiagrGCS(:,2).'],[ElemGCS(:,3).'; SdiagrGCS(:,3).'],'b','Clipping','off');
set(h.hatching,'Color',get(h.stresses,'Color'));
set(h.hatching,'LineStyle',get(h.stresses,'LineStyle'));
set(h.hatching,'Marker',get(h.stresses,'Marker'));

if strcmpi(Values,'on')
    Extrema(find(abs(Extrema)<(max(abs(Extrema))*1e-15)))=0;
    exponent=floor(log10(max(abs(Extrema))))+1;
    if exponent >= 10; fvalues=num2str(Extrema,'%8.5e');
    elseif exponent >= 7;
        for k=1:length(Extrema)
            fvalues(k,:)=sprintf('%7de+3',round(Extrema(k,:)/1000));
        end
    elseif exponent >= 4; fvalues=num2str(Extrema,'%7.0f');
    elseif exponent >= 1; fvalues=num2str(Extrema,'%8.3f'); %fvalues=strrep(fvalues,' 0.000','     0');
    elseif exponent >= -3;
        for k=1:length(Extrema)
            fvalues(k,:)=sprintf('%7de-6\n',round(Extrema(k,:)*1000000));
        end
    else fvalues=num2str(Extrema,'%8.5e');
    end
    h.values=text(ExtremaGCS(:,1),ExtremaGCS(:,2),ExtremaGCS(:,3),fvalues,'HorizontalAlignment','center','VerticalAlignment','bottom');
end

% VIEWPOINT AND AXES FIGURE
if strcmpi(nextplot,'replace')
    if all(Nodes(:,4)==0) && all(SdiagrGCS(:,3)==0 | isnan(SdiagrGCS(:,3)))
        view(0,90);
    elseif all(Nodes(:,3)==0) && all(SdiagrGCS(:,2)==0 | isnan(SdiagrGCS(:,2)))
        view(0,0);
    elseif all(Nodes(:,2)==0) && all(SdiagrGCS(:,1)==0 | isnan(SdiagrGCS(:,1)))
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

if strcmpi(MinMax,'on')
    h.anno=annotation(get(haxis,'Parent'),'textbox',[0.02 0.02 1 1],'LineStyle','none','HorizontalAlignment','left','VerticalAlignment','bottom','String',sprintf('Min: %g    Max: %g',min(Extrema),max(Extrema)),'FontSize',get(0,'defaulttextfontsize'));
    set(get(haxis,'Parent'),'Userdata',h.anno);
end

% RESET NEXTPLOT STATE
set(haxis,'NextPlot',nextplot);

% RETURN OUTPUT ARGUMENTS ONLY IF REQUESTED
if nargout < 1
    clear('StressScal');
end

end