function [Ax,Ay,Az,B,Cx,Cy,Cz,dAxdx,dAydx,dAzdx,dCxdx,dCydx,dCzdx]=disp_rshell(Nodes,Elements,DOF,EltIDDLoad,Sections,Materials,Points,dNodesdx,dDLoadsdx,dSectionsdx)

%DISP_BEAM   Return matrices to compute the displacements of the deformed beams.
%
%   [Ax,Ay,Az,B]=disp_beam(Nodes,Elements,DOF,[],[],[],Points)
%   [Ax,Ay,Az,B]=disp_beam(Nodes,Elements,DOF)
%   returns the matrices to compute the displacements of the deformed beams.
%   The coordinates of the specified points along the deformed beams element are
%   computed using X=Ax*U+B(:,1); Y=Ay*U+B(:,2) and Z=Az*U+B(:,3).
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   DOF        Degrees of freedom  (nDOF * 1)
%   Points     Points in local coordinate system (1 * nPoints)
%   Ax         Matrix to compute the x-coordinates of the deformations
%   Ay         Matrix to compute the y-coordinates of the deformations
%   Az         Matrix to compute the z-coordinates of the deformations
%   B          Matrix which contains the x-, y- and z-coordinates of the
%              undeformed structure
%
%   See also DISP_TRUSS, PLOTDISP.

% Mattias Schevenels
% April 2008

if nargin<8, dNodesdx = []; end
if nargin<9, dDLoadsdx = []; end
if nargin<10, dSectionsdx = []; end
nVar = 0;
if nargout>7 && (~isempty(dNodesdx) || ~isempty(dDLoadsdx) || ~isempty(dSectionsdx))
    nVar = max([size(dNodesdx,3),size(dDLoadsdx,4),size(dSectionsdx,3)]);
end
if nVar>0 && ~isempty(dNodesdx) && nnz(dNodesdx(:,2:end,:))>0
    error('Sensitivities have not been implemented yet.')
end

% PREPROCESSING
DOF=DOF(:);

nDiv=5;
nPoint=4*nDiv+1;
nElt=size(Elements,1);
nDOF=length(DOF);

xi=[linspace(-1,1,nDiv+1),repmat(1,1,nDiv-1),linspace(1,-1,nDiv+1),repmat(-1,1,nDiv)]';
eta=[repmat(-1,1,nDiv),linspace(-1,1,nDiv+1),repmat(1,1,nDiv-1),linspace(1,-1,nDiv+1)]';
zeta=zeros(nPoint,1);

[Nx0,Ny0,Nz0]=nelcs_rshell(xi,eta);

Ax=zeros(nElt*(nPoint+1),nDOF);
Ay=zeros(nElt*(nPoint+1),nDOF);
Az=zeros(nElt*(nPoint+1),nDOF);
B=zeros(nElt*(nPoint+1),3);

for iElt=1:nElt
  eltNode=Nodes(Elements(iElt,5:8),2:4);
  eltNodeNr=Nodes(Elements(iElt,5:8),1);

  Lx=norm(eltNode(2,:)-eltNode(1,:));
  Ly=norm(eltNode(4,:)-eltNode(1,:));

  Tn=trans_beam(eltNode);

  x=[Lx/2*xi,Ly/2*eta,zeta]*Tn'+repmat(mean(eltNode),nPoint,1);

  B((iElt-1)*(nPoint+1)+[1:nPoint],1:3)=x;
  B((iElt-1)*(nPoint+1)+nPoint+1,1:3)=nan;

  Nxe=Nx0;
  Nye=Ny0;
  Nze=Nz0;
  Nze(:,4:6:end)=Nze(:,4:6:end)*Ly/2;
  Nze(:,5:6:end)=Nze(:,5:6:end)*Lx/2;

  Ne=zeros(3*nPoint,24);
  Ne(1:3:end,:)=Nxe;
  Ne(2:3:end,:)=Nye;
  Ne(3:3:end,:)=Nze;

  tmp=repmat({Tn},nPoint,1);
  N=blkdiag(tmp{:})'*Ne*blkdiag(Tn,Tn,Tn,Tn,Tn,Tn,Tn,Tn);

  dofelem=dof_rshell(eltNodeNr);
  C=selectdof(DOF,dofelem);

  Ax((iElt-1)*(nPoint+1)+[1:nPoint],1:nDOF)=N(1:3:end,:)*C;
  Ay((iElt-1)*(nPoint+1)+[1:nPoint],1:nDOF)=N(2:3:end,:)*C;
  Az((iElt-1)*(nPoint+1)+[1:nPoint],1:nDOF)=N(3:3:end,:)*C;
  Ax((iElt-1)*(nPoint+1)+nPoint+1,1:nDOF)=nan;
  Ay((iElt-1)*(nPoint+1)+nPoint+1,1:nDOF)=nan;
  Az((iElt-1)*(nPoint+1)+nPoint+1,1:nDOF)=nan;

  %figure(100);
  %hold('on');
  %plot3(x(:,1),x(:,2),x(:,3),'.-');
  %pause(0.1);


end

if nargout > 4
    Cx = 0;
    Cy = 0; % er wordt geen rekening gehouden met de vervormingen tgv de verdeelde belastingen
    Cz = 0;
    dAxdx = zeros(size(Ax,1),size(Ax,2),nVar);
    dAydx = zeros(size(Ay,1),size(Ay,2),nVar);
    dAzdx = zeros(size(Az,1),size(Az,2),nVar);
    dCxdx = zeros(1,1,nVar);
    dCydx = zeros(1,1,nVar);
    dCzdx = zeros(1,1,nVar);
end

