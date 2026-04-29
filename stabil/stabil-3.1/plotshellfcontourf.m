function h = plotshellfcontourf(ftype,Nodes,Elements,Types,F,varargin)

%PLOTSHELLFCONTOURF   Plot filled contours of shell forces in shell elements.
%
%    plotshellfcontourf(stype,Nodes,Elements,Types,F)
%    plotshellfcontourf(stype,Nodes,Elements,Types,F,DOF,U)
%    plots force contours (output from ELEMSTRESS).
%
%   ftype      'Nx'       Membrane forces 
%              'Ny'     
%              'Nxy'     
%              'Mx'       Bending moments
%              'My'      
%              'Mxy'
%              'Vx'       Shear forces
%              'Vy'
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   F          Forces/moments per unit length 
%              (nElem * 32)              [Nx Ny Nxy Mx My Mxy Vx Vy]
%   DOF        Degrees of freedom  (nDOF * 1)
%   U          Displacements (nDOF * 1)
%   plotshellfcontourf(...,ParamName,ParamValue) sets the value of the
%   specified parameters.  The following parameters can be specified:
%   'GCS'          Plot the GCS. Default: 'on'.
%   'ncolor'       Number of colors in colormap. Default: 10.
%   'Handle'       Plots in the axis with this handle.  Default: current axis.
%
%   See also ELEMSHELLF, PATCH_SHELL8, PATCH_SHELL4.

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
[addminmax,paramlist]=cutparam('minmax','on',paramlist);
[addcbar,paramlist]=cutparam('colorbar','on',paramlist);
paramlist = [{'minmax',addminmax,'colorbar',addcbar},paramlist];


switch lower(ftype)
    case {'nx'}
        ind1 = 0;
    case {'ny'}
        ind1 = 1;
    case {'nxy'}
        ind1 = 2;
    case {'mx'}
        ind1 = 3;
    case {'my'}
        ind1 = 4;
    case {'mxy'}
        ind1 = 5;
    case {'vx'}
        ind1 = 6;
    case {'vy'}
        ind1 = 7;
end

if exist('U','var')
h = plotdispx(Nodes,Elements,Types,DOF,U,F(:,(1:8:25)+ind1),paramlist{:});    
else
h = plotelemx(Nodes,Elements,Types,F(:,(1:8:25)+ind1),paramlist{:});
end

% RETURN OUTPUT ARGUMENTS ONLY IF REQUESTED
if nargout<1, clear('h'); end