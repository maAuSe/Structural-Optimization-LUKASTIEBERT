function [F,dFdx] = loads_solid20(DLoad,Node,dDLoaddx,dNodedx)
%LOADS_SOLID20   Equivalent nodal forces for a solid20 element.
%
%   F = loads_solid20(DLoad,Node)
%   computes the equivalent nodal forces of a distributed load 
%   (in the global coordinate system).
%
%   DLoad      Distributed loads      [n1globalX; n1globalY; n1globalZ; ...]
%              in corner Nodes                       (24 * 1)
%   Node       Node definitions       [x y z] (8 * 3)
%   F          Load vector  (24 * 1)
%
%   See also LOADSLCS_BEAM, ELEMLOADS, LOADS_TRUSS.

% Stijn Francois, Miche Jansen
% 2013

if nargin<3, dDLoaddx = []; end
if nargin<4, dNodedx = []; end
if (nnz(dDLoaddx)+nnz(dNodedx)) > 0
    error('Sensitivities have not been implemented yet.')
end
nVar = 0;
if nargout>1 && (~isempty(dDLoaddx) || ~isempty(dNodedx))
    nVar = max(size(dDLoaddx,3),size(dNodedx,3));
end 

Node = Node(1:20,1:3);
DLoad = DLoad(1:24,:);

% --- Integration points and weights---
nXi1D=2;
xi1D=[5.773502691896258e-001;-5.773502691896258e-001];
H1D =[1;1];

nXi=nXi1D^3;
xi=zeros(nXi,3);
H=zeros(nXi,1);
for iXi=1:nXi1D
  for jXi=1:nXi1D
    for kXi=1:nXi1D
      ind=nXi1D^2*(iXi-1)+nXi1D*(jXi-1)+kXi;
      xi(ind,1)= xi1D(iXi);              % xi1
      xi(ind,2)= xi1D(jXi);              % xi2
      xi(ind,3)= xi1D(kXi);              % xi3
      H(ind)=H1D(iXi)*H1D(jXi)*H1D(kXi);
    end
  end
end

F = zeros(60,size(DLoad,2));

N = zeros(3,60);

for iGauss=1:nXi
    [Ni,dNi_dxi,dNi_deta,dNi_dzeta] = sh_vs20(xi(iGauss,1),xi(iGauss,2),xi(iGauss,3));
    N(1,1:3:22) = Ni;
    N(2,2:3:23) = Ni;
    N(3,3:3:24) = Ni;
    Dloadg = [Ni(1:8,1).'*DLoad(1:3:22,:); Ni(1:8,1).'*DLoad(2:3:23,:); Ni(1:8,1).'*DLoad(3:3:24,:)]; 
    J = Node.'*[dNi_dxi,dNi_deta,dNi_dzeta];
    detJ = det(J);     
    F = F + H(iGauss)*N.'*Dloadg*detJ;
end

if nargout>1
    dFdx = zeros(size(F,1),size(F,2),nVar);  % avoid unassignment error
end
