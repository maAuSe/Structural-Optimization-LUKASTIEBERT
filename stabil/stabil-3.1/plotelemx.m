function h = plotelemx(Nodes,Elements,Types,x,varargin)

%PLOTELEMX  Plot element quantity on elements.
%
%   plotelemx(Nodes,Elements,Types,x) plots element quantity on elements.
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   x          Elements quantity to plot   (nElem * 1)
%   plotelemx(...,ParamName,ParamValue) sets the value of the 
%   specified parameters.  The following parameters can be specified:
%   'GCS'          Plot the GCS. Default: 'on'.
%   'minmax'       Add location of min and max stress. Default: 'off'.
%   'colorbar'     Add colorbar. Default: 'off'.
%   'Undeformed'   Plots the undeformed mesh.  Default: 'k-'.
%   'ncolor'       Number of colors in colormap. Default: 10.
%   'Handle'       Plots in the axis with this handle.  Default: current axis.
%   Additional parameters are redirected to the PATCH function which plots
%   the elements.
%
%   See also ELEMSTRESS.

% Miche Jansen
% 2010

if nargin<5                               
  paramlist={};
elseif nargin>4
  paramlist=varargin;
end

% PARAMETERS
[GCS,paramlist]=cutparam('GCS','on',paramlist);
[addminmax,paramlist]=cutparam('minmax','off',paramlist);
[addcbar,paramlist]=cutparam('colorbar','off',paramlist);
[ncolor,paramlist]=cutparam('ncolor',11,paramlist);
[haxis,paramlist]=cutparam('Handle',gca,paramlist);


[xmax,maxi] = max(max(x,[],2)); 
[xmin,mini] = min(min(x,[],2));

nType = size(Types,1);
pxyz = [];
pind = [];
pvalues = [];

for iType=1:nType
  TypID=Types{iType,1};
  Type=Types{iType,2};
  ind=find(Elements(:,2)==TypID);
  NodeNumbers=Elements(ind,5:end);
  x2 = x(ind,:);
  [pxyz2,pind2,pvalues2]=eval(['patch_' Type '(Nodes,NodeNumbers,x2)']); % aanpassen!!!
  
  if size(pind2,2) < size(pind,2)
    pind2 = [pind2,nan(size(pind2,1),size(pind,2)-size(pind2,2))];
  elseif size(pind,2) < size(pind2,2)
    pind = [pind,nan(size(pind,1),size(pind2,2)-size(pind,2))];
  end
  pind=[pind;pind2+size(pxyz,1)];
  pxyz=[pxyz;pxyz2];
  pvalues=[pvalues;pvalues2];
end

% PREPARE FIGURE
haxis=newplot(haxis);
nextplot=get(haxis,'NextPlot');

% PLOT CONTOUR
% h.elem=patch('Vertices',pxyz,'Faces',pind,'FaceVertexAlphaData',pvalues,'Facealpha','flat','EdgeAlpha','flat',...
%              'FaceVertexCData',pvalues,'Facecolor','flat','EdgeColor','black',paramlist{:});
if size(pind,2) == 2 % line elements
  EdgeC = 'interp';
else
  EdgeC = 'black';
end
h.elem=patch('Vertices',pxyz,'Faces',pind,'FaceVertexCData',pvalues,'Facecolor','interp','EdgeColor',EdgeC,'Clipping','off',paramlist{:});
colormap(jet(ncolor));
set(gcf,'renderer','painters')

% MIN & MAX
if size(x,2) == 1
[loc1,loc2] = ismember(Elements(maxi(1),5:end),Nodes(:,1));
xyzmax = mean(Nodes(loc2(loc1),2:4),1);
[loc1,loc2] = ismember(Elements(mini(1),5:end),Nodes(:,1));
xyzmin = mean(Nodes(loc2(loc1),2:4),1);
else
[xmax,maxj] = max(x(maxi(1),:));
xyzmax = mean(Nodes(Elements(maxi(1),4+maxj(1))==Nodes(:,1),2:4),1);
[xmin,minj] = min(x(mini(1),:));
xyzmin = mean(Nodes(Elements(mini(1),4+minj(1))==Nodes(:,1),2:4),1);
end
set(haxis,'NextPlot','add');
if strcmpi(addminmax,'on')
h.minmax=text([xyzmax(1);xyzmin(1)],[xyzmax(2);xyzmin(2)],[xyzmax(3);xyzmin(3)],{'MAX';'MIN'});
end

% COLORBAR
if strcmpi(addcbar,'on')
v = caxis;
yticks = linspace(v(1),v(2),ncolor+1);
%h.cbar=colorbar('peer',haxis,'YTickMode','manual','location','westoutside','YTickLabelMode','manual','YTick',yticks,'YTickLabel',yticks);
h.cbar=colorbar('peer',haxis,'YTickMode','auto','location','westoutside','YTickLabelMode','auto','YTick',yticks,'YLim',v);
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