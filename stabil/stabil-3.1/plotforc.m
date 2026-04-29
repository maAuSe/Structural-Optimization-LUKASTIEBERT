function [ForcScal,h]=plotforc(ftype,Nodes,Elements,Types,varargin)

%PLOTFORC   Plot element member forces.
%
%   plotforc(ftype,Nodes,Elements,Types,Sections,Materials,Forces,DLoads)
%   plotforc(ftype,Nodes,Elements,Types,Sections,Materials,Forces)
%   plots the element member forces (in accordance to the beam convention).
%
%   ftype      'norm'       Normal force (in the local x-direction)
%   (for truss 'sheary'     Shear force in the local y-direction
%    and beam  'shearz'     Shear force in the local z-direction
%    elements) 'momx'       Torsional moment (around the local x-direction)
%              'momy'       Bending moment around the local y-direction
%              'momz'       Bending moment around the local z-direction
%   ftype      'Nphi'       Normal force (per unit length) in meridional direction
%   (for       'Qphi'       Transverse force (per unit length) in meridional direction
%    shell2    'Mphi'       Bending moment (per unit length) in meridional direction
%    elements) 'Ntheta'     Normal force (per unit length) in circumferential direction
%              'Mtheta'     Bending moment (per unit length) in circumferential direction
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   Sections   Section definitions       [SecID SecProp1 SecProp2 ...]
%   Materials  Material definitions      [MatID MatProp1 MatProp2 ... ]
%   Forces     Element forces in LCS (beam convention) [N Vy Vz T My Mz]
%                                                                   (nElem * 12)
%   DLoads     Distributed loads       [EltID n1globalX n1globalY n1globalZ ...]
%
%   plotforc(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'ForcScal'     Force scaling.  Default: 'auto'.
%   'MinMax'       Display minimum and maximum value.  Default: 'off'.
%   'Values'       Force values.  Default: 'on'.
%   'Undeformed'   Plots the undeformed mesh.  Default: 'k-'.
%   'Handle'       Plots in the axis with this handle.  Default: current axis.
%   Additional parameters are redirected to the PLOT3 function which plots
%   the forces.
%
%   [ForcScal,h] = plotforc(...) returns a struct h with handles to all the
%   objects in the plot.
%
%   See also FDIAGRGCS_BEAM, FDIAGRGCS_TRUSS.

% David Dooms
% October 2008

% Note (MS, April 2020):
% The SHELL2 element needs the Sections and Materials to plot forces along elements;
% the syntax of this function is therefore modified in stabil 3.1:
% - new syntax: plotforc(ftype,Nodes,Elements,Types,Sections,Materials,Forces,DLoads,...)
% - old syntax: plotforc(ftype,Nodes,Elements,Types,Forces,DLoads,...)
% If the old syntax is used, shift the input arguments for backward compatibility
if nargin<7 || isstr(varargin{2}) || isstr(varargin{3})
  Forces = varargin{1};
  varargin = varargin(2:end);
else
  Sections = varargin{1};
  Materials = varargin{2};
  Forces = varargin{3};
  varargin = varargin(4:end);
end

if ~isempty(varargin) && ~isstr(varargin{1})
  DLoads = varargin{1};
  varargin = varargin(2:end);
end
paramlist = varargin;

% PARAMETERS
[ForcScal,paramlist]=cutparam('ForcScal','auto',paramlist);
[MinMax,paramlist]=cutparam('MinMax','off',paramlist);
[Values,paramlist]=cutparam('Values','on',paramlist);
[Undeformed,paramlist]=cutparam('Undeformed','k-',paramlist);
[haxis,paramlist]=cutparam('Handle',gca,paramlist);

ElemGCS=[];
FdiagrGCS=[];
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

    if exist('Sections','var')
      SecID=Elements(iElem,3);
      loc=find(Sections(:,1)==SecID);
      if isempty(loc)
          error('Section %i is not defined.',SecID)
      elseif length(loc)>1
          error('Section %i is multiply defined.',SecID)
      end
      Section=Sections(loc,2:end);
    else
      Section=[];
    end

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
    if ~ exist('DLoads','var')
        DLoads=zeros(1,7,nTimeSteps);
    end
    loc=find(DLoads(:,1)==EltID);
    if isempty(loc)
        DLoad=zeros(6,nTimeSteps);
    elseif length(loc)>=1
        DLoad=permute(DLoads(loc,2:end,:),[2 3 1]);
    end

    [ElemGCS2,FdiagrGCS2,ElemExtGCS2,ExtremaGCS2,Extrema2]=eval(['fdiagrgcs_' Type '(ftype,Force,Node,Section,Material,DLoad,Points)']);

    ElemGCS=[ElemGCS; ElemGCS2; nan(1,3)];
    FdiagrGCS=[FdiagrGCS; FdiagrGCS2; nan(1,3)];
    ElemExtGCS=[ElemExtGCS; ElemExtGCS2];
    ExtremaGCS=[ExtremaGCS; ExtremaGCS2];
    Extrema=[Extrema; Extrema2];

end

% AUTOMATIC FORCE SCALING
if strcmpi(ForcScal,'auto')
    lref=reflength(Nodes);
    ForcScal=3*lref/max(abs(Extrema));
end

FdiagrGCS=ForcScal*FdiagrGCS+ElemGCS;
ExtremaGCS=ForcScal*ExtremaGCS+ElemExtGCS;

% PREPARE FIGURE
nextplot=get(haxis,'NextPlot');

% PLOT MEMBER FORCE
handleforces=plot3(haxis,FdiagrGCS(:,1),FdiagrGCS(:,2),FdiagrGCS(:,3),paramlist{:},'Clipping','off');
set(haxis,'NextPlot','add');

% PLOT UNDEFORMED MESH
if ~ strcmpi(Undeformed,'off')
    h=plotelem(Nodes,Elements,Types,Undeformed,'Numbering','off','Handle',haxis,'GCS','off');
end
h.forces=handleforces;

% PLOT HATCH
h.hatching=plot3([ElemGCS(:,1).'; FdiagrGCS(:,1).'],[ElemGCS(:,2).'; FdiagrGCS(:,2).'],[ElemGCS(:,3).'; FdiagrGCS(:,3).'],'b','Clipping','off');
set(h.hatching,'Color',get(h.forces,'Color'));
set(h.hatching,'LineStyle',get(h.forces,'LineStyle'));
set(h.hatching,'Marker',get(h.forces,'Marker'));

if strcmpi(Values,'on')
    Extrema(find(abs(Extrema)<(max(abs(Extrema))*1e-15)))=0;
    exponent=floor(log10(max(abs(Extrema))))+1;
    if exponent >= 10; fvalues=num2str(Extrema,'%8.5e');
    elseif exponent >= 7;
        for k=1:length(Extrema)
            fvalues{k}=sprintf('%7de+3',round(Extrema(k,:)/1000));
        end
    elseif exponent >= 4; fvalues=num2str(Extrema,'%7.0f');
    elseif exponent >= 0; fvalues=num2str(Extrema,'%8.3f'); %fvalues=strrep(fvalues,' 0.000','     0');
    elseif exponent >= -3;
        for k=1:length(Extrema)
            fvalues{k}=sprintf('%7de-6\n',round(Extrema(k,:)*1000000));
        end
    else fvalues=num2str(Extrema,'%8.5e');
    end
    h.values=text(ExtremaGCS(:,1),ExtremaGCS(:,2),ExtremaGCS(:,3),fvalues,'HorizontalAlignment','center','VerticalAlignment','bottom');
end

% VIEWPOINT AND AXES FIGURE
if strcmpi(nextplot,'replace')
    if all(Nodes(:,4)==0) && all(FdiagrGCS(:,3)==0 | isnan(FdiagrGCS(:,3)))
        view(0,90);
    elseif all(Nodes(:,3)==0) && all(FdiagrGCS(:,2)==0 | isnan(FdiagrGCS(:,2)))
        view(0,0);
    elseif all(Nodes(:,2)==0) && all(FdiagrGCS(:,1)==0 | isnan(FdiagrGCS(:,1)))
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
    h.anno=annotation(get(haxis,'Parent'),'textbox',[0.02 0.02 0.5 0.1],'LineStyle','none','HorizontalAlignment','left','VerticalAlignment','bottom','String',sprintf('Min: %g    Max: %g',min(Extrema),max(Extrema)),'FontSize',get(0,'defaulttextfontsize'));
    set(get(haxis,'Parent'),'Userdata',h.anno);
end

% RESET NEXTPLOT STATE
set(haxis,'NextPlot',nextplot);

% RETURN OUTPUT ARGUMENTS ONLY IF REQUESTED
if nargout < 1
    clear('ForcScal');
end