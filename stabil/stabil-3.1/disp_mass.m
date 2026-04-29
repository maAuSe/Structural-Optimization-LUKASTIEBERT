function [Ax,Ay,Az,B,Cx,Cy,Cz,dAxdx,dAydx,dAzdx,dCxdx,dCydx,dCzdx]=disp_mass(Nodes,Elements,DOF,EltIDDLoad,Sections,Materials,Points,dNodesdx,dDLoadsdx,dSectionsdx)

%DISP_MASS    Matrices to compute the displacements of the deformed mass element
%
%   [Ax,Ay,Az,B,Cx,Cy,Cz]
%            =disp_beam(Nodes,Elements,DOF,EltIDDLoad,Sections,Materials,Points)       
%   [Ax,Ay,Az,B,Cx,Cy,Cz]
%            =disp_beam(Nodes,Elements,DOF,EltIDDLoad,Sections,Materials)
%   [Ax,Ay,Az,B]
%            =disp_beam(Nodes,Elements,DOF,[],Sections,Materials) 
%   [Ax,Ay,Az,B]
%            =disp_beam(Nodes,Elements,DOF)
% 
%   returns the matrices to compute the displacements of the deformed mass.
%   The coordinates of the node of the mass element are 
%   computed using X=Ax*U+B(:,1); Y=Ay*U+B(:,2) and 
%   Z=Az*U+B(:,3). The matrices Cx,Cy and Cz superimpose the 
%   displacements that occur due to the distributed loads if all nodes are fixed.
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   DOF        Degrees of freedom  (nDOF * 1)
%   Ax         Matrix to compute the x-coordinates of the deformations 
%   Ay         Matrix to compute the y-coordinates of the deformations 
%   Az         Matrix to compute the z-coordinates of the deformations 
%   B          Matrix which contains the x-, y- and z-coordinates of the 
%              undeformed structure
%   Cx         Matrix to compute the x-coordinates of the deformations 
%   Cy         Matrix to compute the y-coordinates of the deformations 
%   Cz         Matrix to compute the z-coordinates of the deformations 
%
%   See also DISP_TRUSS, PLOTDISP, DISP_SHELL8.

NodeNum=Elements(:,5).';

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

DOF=DOF(:);
nElem=size(Elements,1);
nNode = 1;

Nx = sparse(1,3);
Ny = sparse(1,3);  
Nz = sparse(1,3);
    
Nx(1,1) = 1;
Ny(1,2) = 1;
Nz(1,3) = 1;

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

X = X(loc);
Y = Y(loc);
Z = Z(loc);

B = [X(:),Y(:),Z(:)];

dofelem = kron(NodeNum(:),ones(3,1))+kron(ones(numel(NodeNum),1),[0.01;0.02;0.03]);

repper = speye(nElem,nElem);
Ax = kron(repper,Nx);
Ay = kron(repper,Ny);
Az = kron(repper,Nz);

Ca=selectdof(DOF,dofelem);
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
    
    
