function [Ax,Ay,Az,B,Cx,Cy,Cz,dAxdx,dAydx,dAzdx,dCxdx,dCydx,dCzdx]=disp_shell2(Nodes,Elements,DOF,EltIDDLoad,Sections,Materials,Points,dNodesdx,dDLoadsdx,dSectionsdx)

%DISP_SHELL2   Return matrices to compute the displacements of deformed SHELL2 elements.
%
%   [Ax,Ay,Az,B,Cx,Cy,Cz]
%            =disp_shell2(Nodes,Elements,DOF,EltIDDLoad,Sections,Materials,Points)
%   [Ax,Ay,Az,B,Cx,Cy,Cz]
%            =disp_shell2(Nodes,Elements,DOF,EltIDDLoad,Sections,Materials)
%   [Ax,Ay,Az,B]
%            =disp_shell2(Nodes,Elements,DOF,[],Sections,Materials)
%   [Ax,Ay,Az,B]
%            =disp_shell2(Nodes,Elements,DOF)
%
%   returns the matrices to compute the displacements of deformed SHELL2 elements.
%   The coordinates of the specified points along the deformed SHELL2 element are
%   computed using X=Ax*U+Cx*DLoad+B(:,1); Y=Ay*U+Cy*DLoad+B(:,2) and
%   Z=Az*U+Cz*DLoad+B(:,3). The matrices Cx,Cy and Cz superimpose the
%   displacements that occur due to the distributed loads if all nodes are fixed.
%   In the current implementation, Cx=Cy=Cz=0, i.e. the local effect of the
%   distributed loads is ignored.
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   DOF        Degrees of freedom  (nDOF * 1)
%   EltIDDLoad Elements with distributed loads [EltID]
%              (use empty array [] when shear deformation is considered 
%                                                         but no DLoads are present)
%   Sections   Section definitions       [SecID SecProp1 SecProp2 ...]
%   Materials  Material definitions      [MatID MatProp1 MatProp2 ... ]
%   Points     Points in the local coordinate system (1 * nPoints)
%   Ax         Matrix to compute the x-coordinates of the deformations
%   Ay         Matrix to compute the y-coordinates of the deformations
%   Az         Matrix to compute the z-coordinates of the deformations
%   B          Matrix which contains the x-, y- and z-coordinates of the
%              undeformed structure
%   Cx         Matrix to compute the x-coordinates of the deformations
%   Cy         Matrix to compute the y-coordinates of the deformations
%   Cz         Matrix to compute the z-coordinates of the deformations

% Mattias Schevenels
% April 2020

% PREPROCESSING
DOF=DOF(:);

if nargin<7
    nPoints=21;
    Points=linspace(0,1,nPoints);
else
    nPoints=length(Points);
end
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

nElem=size(Elements,1);
nDOF=size(DOF,1);
Ax=sparse([],[],[],nElem*(nPoints+1),nDOF,nElem*nPoints*12+nElem);
Ay=sparse([],[],[],nElem*(nPoints+1),nDOF,nElem*nPoints*12+nElem);
Az=sparse([],[],[],nElem*(nPoints+1),nDOF,nElem*nPoints*12+nElem);
B=zeros(nElem*(nPoints+1),3);

NeLCS=nelcs_beam(Points);

for iElem=1:nElem

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

  dofelem=dof_beam(NodeNum);

  C=selectdof(DOF,dofelem);

  % transform displacements from global to local coordinate system
  t = trans_shell2(Node);

  % compute element length
  L=norm(Node(2,:)-Node(1,:));


% Deformations caused by nodal displacements [calculate Ax Ay Az B]
  T=blkdiag(t,t,t,t);
  tempx=NeLCS*diag([1 0 0 0 0 0 1 0 0 0 0 0])*T;
  tempy=NeLCS*diag([0 1 0 0 0 L 0 1 0 0 0 L])*T;
  tempz=NeLCS*diag([0 0 0 0 0 0 0 0 0 0 0 0])*T;

  Ax(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF)=(t(1,1)*tempx+t(2,1)*tempy+t(3,1)*tempz)*C;
  Ay(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF)=(t(1,2)*tempx+t(2,2)*tempy+t(3,2)*tempz)*C;
  Az(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF)=(t(1,3)*tempx+t(2,3)*tempy+t(3,3)*tempz)*C;

  B(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:3)=[NeLCS(:,1), NeLCS(:,7)]*Node(1:2,:);
  B((nPoints+1)*iElem,1)=NaN;
end

% if nargout > 4
%   error('Plotting the effect of distributed loads on SHELL2 element displacements is not supported.');
  % Cx = sparse([],[],[],nElem*(nPoints+1),6*numel(EltIDDLoad));
  % Cy = sparse([],[],[],nElem*(nPoints+1),6*numel(EltIDDLoad));
  % Cz = sparse([],[],[],nElem*(nPoints+1),6*numel(EltIDDLoad));
% end

if nargout>4
    Cx = 0;
    Cy = 0; % er wordt geen rekening gehouden met de vervormingen tgv de verdeelde belastingen
    Cz = 0;
    if nargout>7
        dAxdx = zeros(size(Ax,1),size(Ax,2),nVar);
        dAydx = zeros(size(Ay,1),size(Ay,2),nVar);
        dAzdx = zeros(size(Az,1),size(Az,2),nVar);
        dCxdx = zeros(1,1,nVar);
        dCydx = zeros(1,1,nVar);
        dCzdx = zeros(1,1,nVar);
    end
end

