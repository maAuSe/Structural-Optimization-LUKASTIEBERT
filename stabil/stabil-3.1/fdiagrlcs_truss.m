function [FdiagrLCS,loc,Extrema,dFdiagrLCSdx] = fdiagrlcs_truss(ftype,Forces,~,~,Points,dForcesdx,~,~)

%FDIAGRLCS_TRUSS   Force diagram for a truss element in LCS.
%
%   [FdiagrLCS,loc,Extrema] = FDIAGRLCS_TRUSS(ftype,Forces,DLoadLCS,L,Points)
%   computes the element forces at the specified points. 
%
%   [FdiagrLCS,loc,Extrema,dFdiagrLCSdx] 
%       = FDIAGRLCS_TRUSS(ftype,Forces,DLoadLCS,L,Points,dForcesdx,dDLoadLCSdx,dLdx)
%   additionally computes the derivatives of the element force values with
%   respect to the design variables x.
%
%   ftype           'norm'       Normal force (in the local x-direction)
%   Forces          Element forces in LCS   [N; 0; 0; 0; 0; 0] (12 * 1)
%   Node            Node definitions        [x y z] (3 * 3)
%   Points          Points in the local coordinate system  (1 * nPoints)
%   dForcesdx       Element forces in LCS derivatives      (SIZE(Forces) * nVar)
%   FdiagrLCS       Element forces at the points           (1 * nPoints * nLC)
%   loc             Locations of the extreme values        (nValues * nLC)
%   Extrema         Extreme values                         (nValues * nLC)
%   dFdiagrLCSdx    Element force value derivatives        (1 * nPoints * nLC * nVar)
%
%   See also FDIAGRGCS_TRUSS.

if nargin<6, dForcesdx = []; end

nVar = 0;
if nargout>3 && ~isempty(dForcesdx)
    nVar = size(dForcesdx,3);
end

if nVar==0 || isempty(dForcesdx), dForcesdx = zeros(size(Forces,1),size(Forces,2),nVar); end


Points=Points(:);
nPoints=length(Points);
nLC=size(Forces,2);

dFdiagrLCSdx = zeros(1,nPoints,nLC,nVar);


switch lower(ftype)
    case 'norm'
        FdiagrLCS=permute(Points*Forces(1,:),[3 1 2]);
        for n=1:nVar
            dFdiagrLCSdx(:,:,:,n)=permute(Points*dForcesdx(1,:,n),[3 1 2]);
        end
        loc=repmat([0; 1],1,nLC);
        Extrema=[Forces(1,:); Forces(7,:)];
    case {'sheary','shearz','momx','momy','momz'}
        FdiagrLCS = zeros(1,nPoints,nLC);
        loc=[];
        Extrema=[];
    otherwise
        error('Unknown element force.')
end

