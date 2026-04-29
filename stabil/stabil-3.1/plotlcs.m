function [h,vLCS] = plotlcs(Nodes,Elements,Types,vLCS,varargin)

%PLOTLCS   Plot the local element coordinate systems.
%
%   [h,vLCS] = plotlcs(Nodes,Elements,Types)
%   [h,vLCS] = plotlcs(Nodes,Elements,Types,[],varargin)
%   plotlcs(Nodes,Elements,Types,vLCS,varargin)
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   vLCS       Element coordinate systems (nElem * 9)
%   plotlcs(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'GCS'          Plot the GCS. Default: 'on'.
%   'Undeformed'   Plots the undeformed mesh.  Default: 'k-'.
%   'Handle'       Plots in the axis with this handle.  Default: current
%   axis.
%
%   See also ELEMSTRESS

% Miche Jansen
% 2013

% PARAMETERS
if nargin<5                               
    paramlist={};
elseif nargin>5 && ischar(varargin{1})
    paramlist=varargin;
end


[GCS,paramlist]=cutparam('GCS','on',paramlist);
[Undeformed,paramlist]=cutparam('Undeformed','k:',paramlist);
[haxis,paramlist]=cutparam('Handle',gca,paramlist);


nElem = size(Elements,1);
nType=size(Types,1);
x = zeros(1,nElem);
y = zeros(1,nElem);
z = zeros(1,nElem);

for iType=1:nType
    TypID=Types{iType,1};
    Type=Types{iType,2};
    iElem=find(Elements(:,2)==TypID);
    NodeNumbers=Elements(iElem,5:end);
    [xe,ye,ze]=eval(['coord_' Type '(Nodes,NodeNumbers)']);
    x(iElem) = mean(xe,1);
    y(iElem) = mean(ye,1);
    z(iElem) = mean(ze,1);
end

if nargin < 4 | isempty(vLCS) % LCS not given -> compute LCS
    vLCS = zeros(nElem,9);
    Elements(:,2)=linkandcheck(round(Elements(:,2)),round(cell2mat(Types(:,1))),'Element type');
   [dum,Elements(:,5:end)] = ismember(round(Elements(:,5:end)),round(Nodes(:,1)));
   tfunc = cellfun(@(x) str2func(['trans_' x]),Types(:,2),'UniformOutput',0);
   if size(Types,2) > 2; Options2 = Types(:,3); else Options2 = cell(nType,1); end
   for iElem=1:nElem
       NodeNum=Elements(iElem,5:end);
       Node=NaN(length(NodeNum),3);
       Node(NodeNum~=0,:) = Nodes(NodeNum(NodeNum~=0),2:end);
       t = tfunc{Elements(iElem,2)}(Node,Options2{Elements(iElem,2)}).';
       vLCS(iElem,:) = t(:); 
   end
end

% PREPARE FIGURE
nextplot=get(haxis,'NextPlot');
set(haxis,'NextPlot','add');
handleLCS(1) = quiver3(haxis,x,y,z,vLCS(:,1)',vLCS(:,2)',vLCS(:,3)',0.25,'r');
handleLCS(2) = quiver3(haxis,x,y,z,vLCS(:,4)',vLCS(:,5)',vLCS(:,6)',0.25,'g');
handleLCS(3) = quiver3(haxis,x,y,z,vLCS(:,7)',vLCS(:,8)',vLCS(:,9)',0.25,'b');

% PLOT UNDEFORMED MESH
if ~ strcmpi(Undeformed,'off')
    h=plotelem(Nodes,Elements,Types,Undeformed,'Numbering','off','Handle',haxis,'GCS','off');
end

h.lcs = handleLCS;

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