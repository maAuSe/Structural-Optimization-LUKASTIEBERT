function h = plotstresscontour(stype,Nodes,Elements,Types,S,varargin)

%PLOTSTRESSCONTOUR   Plot stress contour lines in shell elements.
%
%    plotstresscontour(stype,Nodes,Elements,Types,S)
%    plots stress contours (output from ELEMSTRESS/NODALSTRESS).
%
%   stype      'sx'       Normal stress (in the global/local x-direction)
%              'sy'     
%              'sz'     
%              'sxy'      Shear stress 
%              'syz'      
%              'sxz'      
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   S          Element stresses in GCS/LCS  [sxx syy szz sxy syz sxz] 
%                                                                   (nElem * 72)
%   plotstresscontour(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'location'     Location (top,mid,bot). Default: 'top'.
%   'Ncontour'     Number of contours. Default: '10'.
%   'GCS'          Plot the GCS. Default: 'on'.
%   'Undeformed'   Plots the undeformed mesh.  Default: 'k-'.
%   'Handle'       Plots in the axis with this handle.  Default: current axis.
%
%   See also ELEMSTRESS, SCONTOUR_SHELL8, SCONTOUR_SHELL4.

% Miche Jansen
% 2010

if nargin<6                               
    paramlist={};
elseif nargin>5 && ischar(varargin{1})
    paramlist=varargin;
end

% PARAMETERS
[location,paramlist]=cutparam('location','top',paramlist);
[Ncontour,paramlist]=cutparam('Ncontour',10,paramlist);
[GCS,paramlist]=cutparam('GCS','on',paramlist);
[Undeformed,paramlist]=cutparam('Undeformed','k-',paramlist);
[haxis,paramlist]=cutparam('Handle',gca,paramlist);



Scontour=[];
switch lower(stype)
    case {'sx','s3'}
        ind1 = 0;
    case {'sy','s2'}
        ind1 = 1;
    case {'sz','s1'}
        ind1 = 2;
    case {'sxy'}
        ind1 = 3;
    case {'syz'}
        ind1 = 4;
    case {'sxz'}
        ind1 = 5;
end

switch lower(location)
    case {'top'}
        ind2 = 0;
    case {'mid'}
        ind2 = 1;
    case {'bot'};
        ind2 = 2;
end

%[Smax] = max(max(S(:,(1:6:19)+ind1+24*ind2)));
%[Smin] = min(min(S(:,(1:6:19)+ind1+24*ind2)));

[Smax,maxi] = max(S(:,(1:6:19)+ind1+24*ind2)); 
[Smin,mini] = min(S(:,(1:6:19)+ind1+24*ind2));
[Smax,maxj] = max(Smax);
[Smin,minj] = min(Smin);
maxi = maxi(maxj);
mini = mini(minj);

Svalues = linspace(Smin,Smax,Ncontour+2);



nElem=size(Elements,1);

for iElem=1:nElem

    TypID=Elements(iElem,2);
    loc=find(cell2mat(Types(:,1))==TypID);
    if isempty(loc)
        error('Element type %i is not defined.',TypID)
    elseif length(loc)>1
        error('Element type %i is multiply defined.',TypID)
    end
    
    Type=Types{loc,2};
        
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
    
    S2=S(iElem,(1:6:19)+ind1+24*ind2).';
    
Scontour2 =eval(['scontour_' Type '(Node,S2,Svalues)']);
    
    Scontour=[Scontour; Scontour2];
    
end


% PREPARE FIGURE
haxis=newplot(haxis);
nextplot=get(haxis,'NextPlot');

% PLOT CONTOUR
handlestress = [];
teller = 1;
colormap(jet);
while teller < size(Scontour,1)
num = Scontour(teller,2); 
handlestress2=patch('parent',haxis,'xdata',[Scontour((teller+1):(teller+num),1);nan],'ydata',[Scontour((teller+1):(teller+num),2);nan],...
    'zdata',[Scontour((teller+1):(teller+num),3);nan],'facecolor','none','edgecolor','flat','cdata',Scontour(teller,1)+zeros(num+1,1),'LineWidth',2,'Clipping','off');
handlestress = [handlestress;handlestress2];
teller = teller+num+1;

end

% MIN & MAX
loc = find(Nodes(:,1)==Elements(maxi(1),4+maxj(1)));
xyzmax = Nodes(loc(1),2:4);
loc = find(Nodes(:,1)==Elements(mini(1),4+minj(1)));
xyzmin = Nodes(loc(1),2:4);
set(haxis,'NextPlot','add');
h.minmax=text([xyzmax(1);xyzmin(1)],[xyzmax(2);xyzmin(2)],[xyzmax(3);xyzmin(3)],{'MAX';'MIN'});

% COLORBAR
caxis([Smin Smax]);
v = caxis;
yticks = linspace(v(1),v(2),Ncontour+2);
h.cbar=colorbar('peer',haxis,'YTickMode','auto','location','westoutside','YTickLabelMode','auto','YTick',yticks,'YLim',v);


% PLOT UNDEFORMED MESH
if ~ strcmpi(Undeformed,'off')
    %set(haxis,'NextPlot','add');
    h.elem=plotelem(Nodes,Elements,Types,Undeformed,'Numbering','off','Handle',haxis,'GCS','off');
end
h.stress=handlestress;


% VIEWPOINT AND AXES FIGURE
if numel(Scontour) >3
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