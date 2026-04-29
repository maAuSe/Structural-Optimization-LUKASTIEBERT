function [Ax,Ay,Az,B,Cx,Cy,Cz,dAxdx,dAydx,dAzdx,dCxdx,dCydx,dCzdx]=disp_plane10(Nodes,Elements,DOF,EltIDDLoad,Sections,Materials,nPoints,dNodesdx,dDLoadsdx,dSectionsdx)

%DISP_PLANE10   Matrices to compute the displacements of the deformed plane.
%
%   [Ax,Ay,Az,B] = disp_plane10(Nodes,Elements,DOF,U)
%   returns the matrices to compute the displacements of the deformed plane.
%   The coordinates of the nodes of the plane10 element are 
%   computed using X=Ax*U+B(:,1); Y=Ay*U+B(:,2) and 
%   Z=Az*U+B(:,3). 
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   DOF        Degrees of freedom  (nDOF * 1)
%   Ax         Matrix to compute the x-coordinates of the deformations 
%   Ay         Matrix to compute the y-coordinates of the deformations 
%   Az         Matrix to compute the z-coordinates of the deformations 
%   B          Matrix which contains the x-, y- and z-coordinates of the
%              undeformed structure
%
%   See also DISP_TRUSS, PLOTDISP, DISP_SHELL4.

% Miche Jansen, Jef Wambacq
% 2009

DOF=DOF(:);
nElem=size(Elements,1);
nNode = 10;
NodeNum=Elements(:,5:5+(nNode-1)).';

if nargin<8, dNodesdx = []; end
if nargin<9, dDLoadsdx = []; end
if nargin<10, dSectionsdx = []; end
nVar = 0;
if nargout>7 && (~isempty(dNodesdx) || ~isempty(dDLoadsdx) || ~isempty(dSectionsdx))
    nVar = max([size(dNodesdx,3),size(dDLoadsdx,4),size(dSectionsdx,3)]);
end
if nVar>0 && ~isempty(dNodesdx) && nnz(dNodesdx(sort(unique(NodeNum)),2:end,:))>0
    error('Sensitivities have not been implemented yet.')
end

if nargin<7
nPoints=5;
N = zeros(3*(nPoints),nNode);
s = linspace(0,1,nPoints).';
st = [s zeros(nPoints,1);
      1-s s;
      zeros(nPoints,1) flipud(s)]; %nPoints minstens 2 !!!
else
N = zeros(3*(nPoints),nNode);
s = linspace(0,1,nPoints).';
st = [s zeros(nPoints,1);
      1-s s;
      zeros(nPoints,1) flipud(s)];
end

for i = 1:3*(nPoints)
    Ni = sh_t10(st(i,1),st(i,2));
    N(i,:) = Ni; 
end

Nx = sparse(3*nPoints+1,3*nNode);
Ny = sparse(3*nPoints+1,3*nNode);
Nz = sparse(3*nPoints+1,3*nNode);

Nx(1:end-1,1:3:3*nNode-2) = N;
Ny(1:end-1,2:3:3*nNode-1) = N;
Nz(1:end-1,3:3:3*nNode) = N;

% check nodes
if length(unique(Nodes(:,1)))~=length(Nodes(:,1))
    for iNode = 1:length(Nodes(:,1))
        loc = find(abs(Nodes(iNode,1)-Nodes(:,1))<1e-4);
        if length(loc)>1
        error('Node %i is multiply defined.',Nodes(iNode,1))
        end
    end
end

X = Nodes(:,2);
Y = Nodes(:,3);
Z = Nodes(:,4);

[exis,loc] = ismember(round(NodeNum),round(Nodes(:,1)));

if ~all(exis(:))
   unknowns = NodeNum(~exis);
   error('Node %i is not defined.',unknowns(1));
end

X = [N*X(loc);nan(1,size(loc,2))];
Y = [N*Y(loc);nan(1,size(loc,2))];
Z = [N*Z(loc);nan(1,size(loc,2))];

B = [X(:),Y(:),Z(:)];

% select translation-dofs
dofelem = kron(NodeNum(:),ones(3,1))+kron(ones(numel(NodeNum),1),[0.01;0.02;0.03]);

Ca=selectdof(DOF,dofelem);

repper = speye(nElem,nElem);
Ax = kron(repper,Nx);
Ay = kron(repper,Ny);
Az = kron(repper,Nz);
   
Ax = Ax*Ca;
Ay = Ay*Ca;
Az = Az*Ca;

if nargout > 4
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
end
