function h = plotstresscontourf(stype,Nodes,Elements,Types,S,varargin)

%PLOTSTRESSCONTOURF  Plot filled contours of stresses.
%
%    plotstresscontourf(stype,Nodes,Elements,Types,S)
%    plotstresscontourf(stype,Nodes,Elements,Types,S,DOF,U)
%    plots stress contours (output from ELEMSTRESS).
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
%                                                              (nElem * 72)
%   plotstresscontourf(...,ParamName,ParamValue) sets the value of the 
%   specified parameters.  The following parameters can be specified:
%   'location'     Location (top,mid,bot). Default: 'top'.
%   'GCS'          Plot the GCS. Default: 'on'.
%   'minmax'       Add location of min and max stress. Default: 'on'.
%   'colorbar'     Add colorbar. Default: 'on'.
%   'ncolor'       Number of colors in colormap. Default: 10.
%   'Handle'       Plots in the axis with this handle.  Default: current axis.
%
%   See also ELEMSTRESS.

% Miche Jansen
% 2010


if nargin<6                               
    paramlist={};
elseif nargin>5 && ischar(varargin{1})
    paramlist=varargin;
elseif nargin>6 && ischar(varargin{2})
    paramlist=varargin(2:end);
    warning('The sixth input argument is not used! Check syntax.')
elseif nargin==7
    DOF=varargin{1};
    U=varargin{2};
    paramlist = {};
else
    DOF=varargin{1};
    U=varargin{2};
    paramlist = varargin(3:end);
end

% PARAMETERS
[location,paramlist]=cutparam('location','top',paramlist);
[addminmax,paramlist]=cutparam('minmax','on',paramlist);
[addcbar,paramlist]=cutparam('colorbar','on',paramlist);
paramlist = [{'minmax',addminmax,'colorbar',addcbar},paramlist];

S2 = zeros(size(Elements,1),1);
for iType = 1:size(Types,1)
    TypID=Types{iType,1};
    Type=Types{iType,2};
    loc=find(Elements(:,2)==TypID);
    es = eltdef(Type); 
    if isstruct(es.stress.(stype))
    S2(loc,1:length(es.stress.(stype).(location))) = S(loc,es.stress.(stype).(location));
    else
    S2(loc,1:length(es.stress.(stype))) = S(loc,es.stress.(stype));        
    end
end

if exist('U','var')
h = plotdispx(Nodes,Elements,Types,DOF,U,S2,paramlist{:});    
else
h = plotelemx(Nodes,Elements,Types,S2,paramlist{:});
end

% RETURN OUTPUT ARGUMENTS ONLY IF REQUESTED
if nargout<1, clear('h'); end