function plotelemsec(Nodes,Elements,Types,Sections,varargin)

%PLOTELEMSEC   Plot elements with line thicknesses proportional to section areas.
%
%   PLOTELEMSEC(Nodes,Elements,Types,Sections)
%   plots the elements.
%
%   Nodes      Node definitions          [NodID x y z]
%   Elements   Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   Sections   Section definitions       [SecID SecProp1 SecProp2 ...]
%
%   PLOTELEMSEC(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'Numbering'        Plots the element numbers.  Default: 'on'.
%   'GCS'              Plots the global coordinate system.  Default: 'on'.
%   'Handle'           Plots in the axis with this handle.  Default: current axis.
%   'LineWidthScale'   Set line width scale.  Default: 8/MAX(A).
%   Additional parameters are redirected to the PLOT3 function which plots
%   the elements.
%
%   See also PLOTELEM.

% Mattias Schevenels
% February 2017

SecIDs = unique(Elements(:,3));
A = Sections(:,2);

paramlist = varargin;
[LineWidthScale,paramlist]=cutparam('LineWidthScale',8/max(A),paramlist);

for k = 1:length(SecIDs)
  plotelem(Nodes,Elements(Elements(:,3)==SecIDs(k),:),Types,paramlist{:},'LineWidth',A(k)*LineWidthScale)
  if k==1,
    oldNextPlot = get(gca,'NextPlot');
  end
  hold('on');
end
set(gca,'NextPlot',oldNextPlot);


function [value,paramlist]=cutparam(name,default,paramlist);
% CUT PARAMETER FROM LIST
value=default;
for iarg=length(paramlist)-1:-1:1
    if strcmpi(name,paramlist{iarg})
        value=paramlist{iarg+1};
        paramlist=paramlist([1:iarg-1 iarg+2:end]);
        break
    end
end