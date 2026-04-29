function [DispScal,h,Ax,Ay,Az,B,Cx,Cy,Cz]=plotdisp(Nodes,Elements,Types,DOF,U,varargin) 

%PLOTDISP   Plot the displacements.
% 
%            plotdisp(Nodes,Elements,Types,DOF,U,DLoads,Sections,Materials)
%            plotdisp(Nodes,Elements,Types,DOF,U,[],Sections,Materials)
%            plotdisp(Nodes,Elements,Types,DOF,U)
%   DispScal=plotdisp(Nodes,Elements,Types,DOF,U,DLoads,Sections,Materials)
%   plots the displacements. If DLoads, Sections and Materials are supplied, the 
%   displacements that occur due to the distributed loads if all nodes are fixed,
%   are superimposed.
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   DOF        Degrees of freedom  (nDOF * 1)
%   U          Displacements (nDOF * 1)
%   DLoads     Distributed loads       [EltID n1globalX n1globalY n1globalZ ...]
%              (use empty array [] when shear deformation (in beam element) 
%               is considered but no DLoads are present)
%   Sections   Section definitions     [SecID SecProp1 SecProp2 ...]
%   Materials  Material definitions    [MatID MatProp1 MatProp2 ... ]
%   DispScal   Displacement scaling
%
%   plotdisp(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'DispScal'     Displacement scaling.  Default: 'auto'.
%   'DispMax'      Mention maximal displacement.  Default: 'on'.
%   'Undeformed'   Plots the undeformed mesh.  Default: 'k:'.
%   'Handle'       Plots in the axis with this handle.  Default: current axis.
%   Additional parameters are redirected to the PLOT3 function which plots 
%   the deformations.
%
%   [DispScal,h] = plotdisp(...) returns a struct h with handles to all the 
%   objects in the plot.
%
%   See also DISP_TRUSS, DISP_BEAM, PLOTELEM.

% David Dooms
% September 2008

if nargin<6                             % plotdisp(Nodes,Elements,Types,DOF,U)
    paramlist={};
elseif nargin>5 && ischar(varargin{1})
    paramlist=varargin;
elseif nargin>6 && ischar(varargin{2})
    paramlist=varargin(2:end);
    warning('The sixth input argument is not used! Check syntax.') 
elseif nargin>7 && ischar(varargin{3})
    paramlist=varargin(3:end);
    warning('The seventh input argument is not used! Check syntax.')
elseif nargin==8            % plotdisp(Nodes,Elements,Types,DOF,U,DLoads,Sections,Materials)
    DLoads=varargin{1};
    Sections=varargin{2};
    Materials=varargin{3};
    paramlist={}; 
else
    DLoads=varargin{1};
    Sections=varargin{2};
    Materials=varargin{3};
    paramlist=varargin(4:end); 
end

% PARAMETERS
[DispScal,paramlist]=cutparam('DispScal','auto',paramlist);
[DispMax,paramlist]=cutparam('DispMax','on',paramlist);
[Undeformed,paramlist]=cutparam('Undeformed','k:',paramlist);
[haxis,paramlist]=cutparam('Handle',gca,paramlist);

nType=size(Types,1);

Ax=[];
Ay=[];
Az=[];
B=[];

if exist('DLoads','var') && ~isempty(DLoads)
    Cx=[];
    Cy=[];
    Cz=[];
end

for iType=1:nType
    loc=find(Elements(:,2)==cell2mat(Types(iType,1)));
    if isempty(loc)
        
    else
        ElemType=Elements(loc,:);
        Type=Types{iType,2};
        
        if exist('DLoads','var')
            if isempty(DLoads)
            [Axadd,Ayadd,Azadd,Badd]=eval(['disp_' Type '(Nodes,ElemType,DOF,[],Sections,Materials)']);
            else
            [Axadd,Ayadd,Azadd,Badd,Cxadd,Cyadd,Czadd]=eval(['disp_' Type '(Nodes,ElemType,DOF,DLoads,Sections,Materials)']);
            Cx=[Cx; Cxadd];
            Cy=[Cy; Cyadd];
            Cz=[Cz; Czadd];
            end
        else
            [Axadd,Ayadd,Azadd,Badd]=eval(['disp_' Type '(Nodes,ElemType,DOF)']);
        end 
        Ax=[Ax; Axadd];
        Ay=[Ay; Ayadd];
        Az=[Az; Azadd];
        B=[B; Badd];
    end
end

nTSteps=size(U,2);

if exist('DLoads','var') && ~isempty(DLoads)
    DLoads=permute(DLoads(:,2:7,:),[2 1 3]);
    DLoad=reshape(DLoads,numel(DLoads(:,:,1)),size(DLoads,3));
    X=Ax*U+Cx*DLoad;
    Y=Ay*U+Cy*DLoad;
    Z=Az*U+Cz*DLoad;
else
    X=Ax*U;
    Y=Ay*U;
    Z=Az*U;
end

% CHECK SOLUTION
lref=reflength(Nodes);
DispAmpl=sqrt(X.^2+Y.^2+Z.^2);
MaxDisp=max(max(DispAmpl));
%%if MaxDisp>lref*1e9
%%    warning('The displacements are very large. Check for an insufficiently constrained model. Press any key to continue.');
%%    pause
%%elseif any(any(isnan(DispAmpl))) 
%%    warning('Not all of the displacements are real. Check for an insufficiently constrained model. Press any key to continue.');
%%    pause
%%end

% AUTOMATIC DISPLACEMENT SCALING
if strcmpi(DispScal,'auto')
    DispScal=3*lref/MaxDisp;
end

X=DispScal*X+repmat(B(:,1),1,nTSteps);
Y=DispScal*Y+repmat(B(:,2),1,nTSteps);
Z=DispScal*Z+repmat(B(:,3),1,nTSteps);

% PREPARE FIGURE
nextplot=get(haxis,'NextPlot');

% PLOT DEFORMED MESH
handledisp=plot3(haxis,X,Y,Z,paramlist{:},'Clipping','off');

% PLOT UNDEFORMED MESH
if ~ strcmpi(Undeformed,'off')
    set(haxis,'NextPlot','add');
    h=plotelem(Nodes,Elements,Types,Undeformed,'Numbering','off','Handle',haxis,'GCS','off');
end
h.displacements=handledisp;

% VIEWPOINT AND AXES FIGURE
if strcmpi(nextplot,'replace')
    if all(Nodes(:,4)==0) && all(all((selectdof(DOF,[0.03; 0.04; 0.05])*U)==0))
        view(0,90);
    elseif all(Nodes(:,3)==0) && all(all((selectdof(DOF,[0.02; 0.04; 0.06])*U)==0))
        view(0,0);
    elseif all(Nodes(:,2)==0) && all(all((selectdof(DOF,[0.01; 0.05; 0.06])*U)==0))
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
    %% axis fill  % octave compatibility
    axis off
    set(gcf,'PaperPositionMode','auto')
end

if strcmpi(DispMax,'on')
    h.anno=annotation(get(haxis,'Parent'),'textbox',[0.02 0.02 1 1],'LineStyle','none','HorizontalAlignment','left','VerticalAlignment','bottom','String',sprintf('Maximal displacement: %g',full(MaxDisp)),'FontSize',get(0,'defaulttextfontsize'));
    set(get(haxis,'Parent'),'Userdata',h.anno);
end

% RESET NEXTPLOT STATE
set(haxis,'NextPlot',nextplot);

% RETURN OUTPUT ARGUMENTS ONLY IF REQUESTED
if nargout < 1
    clear('DispScal');
end
end