function [ForcesLCS,ForcesGCS,dForcesLCSdx,dForcesGCSdx]=forces_shell2(Node,Section,Material,UeGCS,DLoad,TLoad,Options,dNodedx,dSectiondx,dUeGCSdx,dDLoaddx)

%FORCES_BEAM   Compute the element forces for a SHELL2 element.
%
%   [ForcesLCS,ForcesGCS]=forces_shell2(Node,Section,Material,UeGCS,DLoad,TLoad,Options)
%   [ForcesLCS,ForcesGCS]=forces_shell2(Node,Section,Material,UeGCS,DLoad)
%   [ForcesLCS,ForcesGCS]=forces_shell2(Node,Section,Material,UeGCS)
%   computes the element forces for the SHELL2 element in the local and the
%   global coordinate system (algebraic convention).
%
%   Node       Node definitions           [x y z] (3 * 3)
%   Section    Section definition         [h]
%   Material   Material definition        [E nu]
%   UeGCS      Displacements (12 * 1)
%   DLoad      Distributed loads       [n1globalX; n1globalY; n1globalZ; ...]
%                                                                        (6 * 1)
%   TLoad      TLoad                      [dTm; dTy; dTz] (3 * 1)
%   Options    Element options            {Option1 Option2 ...}
%   ForcesLCS  Element forces in the LCS  (12 * 1)
%   ForcesGCS  Element forces in the GCS  (12 * 1)

% Mattias Schevenels
% April 2020

if nargin<5, DLoad=zeros(6,size(UeGCS,2)); end
if nargin<6, TLoad = []; end
if nargin<8, dNodedx = []; end
if nargin<9, dSectiondx = []; end
if nargin<10, dUeGCSdx = []; end
if nargin<11, dDLoaddx = []; end
if (nnz(dNodedx)+nnz(dSectiondx)+nnz(dUeGCSdx)+nnz(dDLoaddx)) > 0
    error('Sensitivities have not been implemented yet.')
end
nVar = 0;
if nargout>2 && (~isempty(dNodedx) || ~isempty(dSectiondx) || ~isempty(dUeGCSdx) || ~isempty(dDLoaddx))
    nVar = max([size(dNodedx,3),size(dSectiondx,3),size(dUeGCSdx,3),size(dDLoaddx,3)]);
end 

dForcesLCSdx=zeros(size(UeGCS,1),size(UeGCS,2),nVar);
dForcesGCSdx=zeros(size(UeGCS,1),size(UeGCS,2),nVar);

% TRANSFORM DISPLACEMENTS FROM GLOBAL TO LOCAL COORDINATE SYSTEM
t=trans_shell2(Node);
T=blkdiag(t,t);
UeLCS=T*UeGCS;

% TRANSFORM DISTRIBUTED LOADS FROM GLOBAL TO LOCAL COORDINATE SYSTEM
if size(DLoad,1)==12 % distributed forces and moments
  T=blkdiag(t,t,t,t);
  DLoadLCS=T*DLoad(1:12,:);
else % only distributed forces
  T=blkdiag(t,t);
  DLoadLCS=T*DLoad(1:6,:);
end

% NO TEMPERATURE LOADS
if nargin>5 && any(TLoad(:)~=0)
  error('Determining the effect of temperature loads on internal SHELL2 element forces is not supported.');
end

% ELEMENT LENGTH, LOCATION, SLOPE
L = norm(Node(2,:)-Node(1,:));
r1 = Node(1,1);
phi = atan2(Node(2,2)-Node(1,2),Node(2,1)-Node(1,1));

E=Material(1);
nu=Material(2);
h=Section(1);

% COMPUTE THE ELEMENT FORCES IN THE LOCAL COORDINATE SYSTEM
KeLCS = kelcs_shell2(r1,phi,L,h,E,nu);
ForcesLCS = KeLCS*UeLCS;

% SUBTRACT EQUIVALENT NODAL FORCES (DISTRIBUTED LOADS) FROM ELEMENT FORCES
ForcesLCS = ForcesLCS-loadslcs_shell2(DLoadLCS,L);

% TRANSFORM THE ELEMENT FORCES FROM LOCAL TO GLOBAL COORDINATE SYSTEM
T=blkdiag(t,t);
ForcesGCS=T.'*ForcesLCS;

% ADD (ZERO) OUT-OF-PLANE COMPONENTS
ForcesLCS = [ForcesLCS(1:2); 0; 0; 0; ForcesLCS(3:5); 0; 0; 0; ForcesLCS(6)];
ForcesGCS = [ForcesGCS(1:2); 0; 0; 0; ForcesGCS(3:5); 0; 0; 0; ForcesGCS(6)];
