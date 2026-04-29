function [Ax,Ay,Az,B,Cx,Cy,Cz,dAxdx,dAydx,dAzdx,dCxdx,dCydx,dCzdx]=disp_truss(Nodes,Elements,DOF,DLoads,Sections,Materials,Points,dNodesdx,dDLoadsdx,dSectionsdx)

%DISP_TRUSS  Return matrices to compute the displacements of the deformed trusses.
%
%   [Ax,Ay,Az,B]=disp_truss(Nodes,Elements,DOF,[],[],[],Points)
%   [Ax,Ay,Az,B]=disp_truss(Nodes,Elements,DOF)
%       returns the matrices to compute the displacements of the deformed
%       trusses. The coordinates of the specified points along the deformed
%       truss elements are computed using X=Ax*U+B(:,1); Y=Ay*U+B(:,2) 
%       and Z=Az*U+B(:,3).
%
%   [Ax,Ay,Az,B,Cx,Cy,Cz,dAxdx,dAydx,dAzdx,dCxdx,dCydx,dCzdx] 
%           = DISP_TRUSS(Nodes,Elements,DOF,EltIDDLoad,Sections,Materials,Points,dNodesdx,
%                                                                       dDLoadsdx(,dSectionsdx))
%       additionally computes the derivatives of the displacements with
%       respect to the design variables x.
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   DOF        Degrees of freedom  (nDOF * 1)
%   Points     Points in the local coordinate system  (1 * nPoints)
%   dNodesdx   Node definitions derivatives   (SIZE(Node) * nVar)
%   dDLoads    Distributed loads derivatives  (SIZE(DLoad) * nVar)
%   Ax         Matrix to compute the x-coordinates of the deformations 
%   Ay         Matrix to compute the y-coordinates of the deformations 
%   Az         Matrix to compute the z-coordinates of the deformations 
%   B          Matrix which contains the x-, y- and z-coordinates of the 
%              undeformed structure
%   dAxdx, dAydx, dAzdx, dCxdx, dCydx, dCzdx    
%       Derivatives of the matrices to compute the coordinates of the
%       interpolation points after deformation
%
%   See also DISP_BEAM, PLOTDISP.

% David Dooms, Wouter Dillen
% March 2008, July 2017

if nargin<4, DLoads = []; end
if nargin<7, Points = []; end
if nargin<8, dNodesdx = []; end
if nargin<9, dDLoadsdx = []; end
if nargin<10, dSectionsdx = []; end

nVar = 0;
if nargout>7 && (~isempty(dNodesdx) || ~isempty(dDLoadsdx) || ~isempty(dSectionsdx))
    nVar = max([size(dNodesdx,3),size(dDLoadsdx,4),size(dSectionsdx,3)]);
end

if isempty(Points)
    nPoints=2;
    Points = linspace(0,1,nPoints);
else
    Points = Points(:).';
    nPoints=length(Points);
end


nElem=size(Elements,1);
nDLoad=size(DLoads,1);
DOF=DOF(:);
nDOF=size(DOF,1);

Ax=sparse([],[],[],nElem*(nPoints+1),nDOF,nElem*nPoints*12+nElem);
Ay=sparse([],[],[],nElem*(nPoints+1),nDOF,nElem*nPoints*12+nElem);
Az=sparse([],[],[],nElem*(nPoints+1),nDOF,nElem*nPoints*12+nElem);
B=zeros(nElem*(nPoints+1),3);
%
dAxdx=zeros(nElem*(nPoints+1),nDOF,nVar);
dAydx=zeros(nElem*(nPoints+1),nDOF,nVar);
dAzdx=zeros(nElem*(nPoints+1),nDOF,nVar);

Cx=sparse([],[],[],nElem*(nPoints+1),nDLoad*6,0);
Cy=sparse([],[],[],nElem*(nPoints+1),nDLoad*6,0);
Cz=sparse([],[],[],nElem*(nPoints+1),nDLoad*6,0);
%
dCxdx=zeros(nElem*(nPoints+1),nDLoad*6,nVar);
dCydx=zeros(nElem*(nPoints+1),nDLoad*6,nVar);
dCzdx=zeros(nElem*(nPoints+1),nDLoad*6,nVar);


% Interpolation functions
A=[-1  1;
   -1  1;
   -1  1;
    1  0;
    1  0;
    1  0];

NeLCS=zeros(length(Points),2);
for k=1:6
    NeLCS(:,k)=polyval(A(k,:),Points);
end


% Loop over the elements
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
    
    dofelem=dof_truss(NodeNum);
    C=selectdof(DOF,dofelem);
    
    Ax(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF)=NeLCS*diag([1 0 0 1 0 0])*C;
    Ay(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF)=NeLCS*diag([0 1 0 0 1 0])*C;
    Az(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF)=NeLCS*diag([0 0 1 0 0 1])*C;
    
    B(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:3)=[NeLCS(:,1), NeLCS(:,4)]*Node(1:2,:);
    B((nPoints+1)*iElem,1)=NaN;

end

