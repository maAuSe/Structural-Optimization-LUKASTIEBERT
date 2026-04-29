function [VectScal,h] = plotprincstress(Nodes,Elements,Types,Spr,Vpr,varargin)

%PLOTPRINCSTRESS   Plot the principal stresses in shell elements.
%
%   plotprincstress(Nodes,Elements,Types,Spr,Vpr)
%   plots the principal stresses in shell elements with a vector plot.
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   Spr        Principal stresses (nElem * 72) [s3 s2 s1 0 0 0 ...]
%   Vpr        Principal dir. matrix (output from principalstress)
%              {nElem * 12}
%   plotprincdir(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'location'     Location (top,mid,bot). Default: 'top'.
%   'GCS'          Plot the GCS. Default: 'on'.
%   'VectScal'     Vector scaling.  Default: 'auto'.
%   'Undeformed'   Plots the undeformed mesh.  Default: 'k-'.
%   'Handle'       Plots in the axis with this handle.  Default: current
%   axis.
%
%   See also PRINCIPALSTRESS.

% PARAMETERS
if nargin<6                               
    paramlist={};
elseif nargin>6 && ischar(varargin{1})
    paramlist=varargin;
end

[location,paramlist]=cutparam('location','top',paramlist);
[GCS,paramlist]=cutparam('GCS','on',paramlist);
[Undeformed,paramlist]=cutparam('Undeformed','k:',paramlist);
[VectScal,paramlist]=cutparam('VectScal','auto',paramlist);
[haxis,paramlist]=cutparam('Handle',gca,paramlist);

switch lower(location)
    case {'top'}
        iloc = 0;
    case {'mid'}
        iloc = 1;
    case {'bot'};
        iloc = 2;
end

nElem = size(Elements,1);
x = [];
y = [];
z = [];
v3 = [];
v2 = [];
v1 = [];
s1 = [];
s2 = [];
s3 = [];

for iElem=1:nElem
    
    TypID=Elements(iElem,2);
    loc=find(cell2mat(Types(:,1))==TypID);
    if isempty(loc)
        error('Element type %i is not defined.',TypID)
    elseif length(loc)>1
        error('Element type %i is multiply defined.',TypID)
    end
    
    Type=Types{loc,2};
    
if any(strcmpi(Type,{'shell8','shell4','shell6'}))
    switch lower(Type)
        case {'shell8','shell4'} % Enkel hoekknopen worden beschouwd
        nNode = 4;
        case 'shell6'
        nNode = 3;    
    end 
    NodeNum=Elements(iElem,4+(1:nNode));
    
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
    

% elementen moeten vier knopen hebben voor deze functie !!
x2 = Node(1:nNode,1)';
y2 = Node(1:nNode,2)';
z2 = Node(1:nNode,3)';


va = [Vpr{iElem,4*(iloc)+(1:nNode)}];

x = [x,x2];
y = [y,y2];
z = [z,z2];


s1a = Spr(iElem,3+(0:6:(6*(nNode-1)))+24*iloc);
s2a = Spr(iElem,2+(0:6:(6*(nNode-1)))+24*iloc);
s3a = Spr(iElem,1+(0:6:(6*(nNode-1)))+24*iloc);

s1 = [s1,s1a];
s2 = [s2,s2a];
s3 = [s3,s3a];

v3 = [v3,va(:,1+(0:3:3*(nNode-1)))];
v2 = [v2,va(:,2+(0:3:3*(nNode-1)))];
v1 = [v1,va(:,3+(0:3:3*(nNode-1)))];
end

end

v1 = v1.*[s1;s1;s1];
v2 = v2.*[s2;s2;s2];
v3 = v3.*[s3;s3;s3];

if strcmpi(VectScal,'auto')
% scaling based on the procedure in quiver3
xPoints = length(unique(x));
yPoints = length(unique(y));
zPoints = length(unique(z));

delx = (min(x)-max(x))/xPoints;
dely = (min(y)-max(y))/yPoints;
delz = (min(z)-max(z))/zPoints;
del = sqrt(delx^2+dely^2+delz^2);
VectScal = sqrt(([v1(1,:),v2(1,:),v3(1,:)]).^2+([v1(2,:),v2(2,:),v3(2,:)]).^2+([v1(3,:),v2(3,:),v3(3,:)]).^2)/del;
VectScal = max(VectScal);
VectScal = 1./VectScal;

v1=0.45*v1*VectScal;
v2= 0.45*v2*VectScal;
v3= 0.45*v3*VectScal;
else v1=v1*VectScal;v2=v2*VectScal;v3=v3*VectScal; 
end


V1=[v1,-v1];V2=[v2,-v2];V3=[v3,-v3];S1=[s1,s1];S2=[s2,s2];S3=[s3,s3];

% PREPARE FIGURE
nextplot=get(haxis,'NextPlot');
set(haxis,'NextPlot','add');

X=[x,x];Y=[y,y];Z=[z,z];
X(S3<0)= X(S3<0) + V3(1,S3<0);
Y(S3<0)= Y(S3<0) + V3(2,S3<0);
Z(S3<0)= Z(S3<0) + V3(3,S3<0);
V3(:,S3<0) = -V3(:,S3<0);
handleprinccs(1) = quiver3(haxis,X,Y,Z,V3(1,:),V3(2,:),V3(3,:),0,'b');

X=[x,x];Y=[y,y];Z=[z,z];
X(S2<0)= X(S2<0) + V2(1,S2<0);
Y(S2<0)= Y(S2<0) + V2(2,S2<0);
Z(S2<0)= Z(S2<0) + V2(3,S2<0);
V2(:,S2<0) = -V2(:,S2<0);
handleprinccs(2) = quiver3(haxis,X,Y,Z,V2(1,:),V2(2,:),V2(3,:),0,'b');

X=[x,x];Y=[y,y];Z=[z,z];
X(S1<0)= X(S1<0) + V1(1,S1<0);
Y(S1<0)= Y(S1<0) + V1(2,S1<0);
Z(S1<0)= Z(S1<0) + V1(3,S1<0);
V1(:,S1<0) = -V1(:,S1<0);
handleprinccs(3) = quiver3(haxis,X,Y,Z,V1(1,:),V1(2,:),V1(3,:),0,'b');

% PLOT UNDEFORMED MESH
if ~ strcmpi(Undeformed,'off')
    h=plotelem(Nodes,Elements,Types,Undeformed,'Numbering','off','Handle',haxis,'GCS','off');
end

h.princcs = handleprinccs;

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