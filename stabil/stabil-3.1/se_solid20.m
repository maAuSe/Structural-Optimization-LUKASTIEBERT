function [SeGCS,SeLCS,vLCS] = se_solid20(Node,Section,Material,UeGCS,Options,gcs)

%SE_SOLID20   Compute the element stresses for a solid20 element.
%
%   [SeGCS,SeLCS,vLCS] = se_solid20(Node,Section,Material,UeGCS,Options,gcs)
%   [SeGCS,SeLCS]      = se_solid20(Node,Section,Material,UeGCS,Options,gcs)
%    SeGCS             = se_solid20(Node,Section,Material,UeGCS,Options,gcs)
%   computes the element stresses in the global and the
%   local coordinate system for the solid20 element.
%
%   Node       Node definitions           [x y z] (20 * 3)
%              Nodes should have the following order:
%                   8---15----7
%                  /|        /|
%                16 |      14 |
%                / 20      / 19 
%               /   |     /   |
%              5---13----6    |
%              |    4--11|----3
%              |   /     |   /
%             17 12     18  10
%              | /       | /
%              |/        |/
%              1----9----2
%

%   Section    Section definition         []
%   Material   Material definition        [E nu rho]
%   UeGCS      Displacements (48 * nTimeSteps)
%   Options    Element options struct. Fields: []
%   GCS        Global coordinate system in which stresses are returned
%              'cart'|'cyl'|'sph'
%   SeGCS      Element stresses in GCS in nodes (48 * nTimeSteps)
%              48 = 6 stress comp. * 8 corner nodes
%                                        [sxx syy szz sxy syz sxz]
%   SeGCS      Element stresses in GCS in nodes (48 * nTimeSteps)
%              48 = 6 stress comp. * 8 corner nodes
%                                        [sxx syy szz sxy syz sxz]
%   vLCS       Unit vectors of LCS (3 * 3)
%
%   See also ELEMSTRESS, SE_SHELL4.

% Stijn Francois, Miche Jansen
% 2013

Node = Node(1:20,:);

% INITIALIZATION
Be=zeros(6,60); % strain matrix
SeGCS=zeros(48,size(UeGCS,2));

C=cmat_isotropic('3d',Section,Material);

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

% stress in Gauss points
SeGCSg=zeros(6*nXi,size(UeGCS,2));
iudof = (1:3:58);
ivdof = 1+iudof;
iwdof = 2+iudof;

for iGauss=1:nXi
    
    [Ni,dNi_dxi,dNi_deta,dNi_dzeta] = sh_vs20(xi(iGauss,1),xi(iGauss,2),xi(iGauss,3));
    
    J = Node.'*[dNi_dxi,dNi_deta,dNi_dzeta];
    dNi_dX = J\[dNi_dxi,dNi_deta,dNi_dzeta].';
    
    Be(1,iudof) = dNi_dX(1,:);
    Be(2,ivdof) = dNi_dX(2,:);
    Be(3,iwdof) = dNi_dX(3,:);
    Be(4,iudof) = dNi_dX(2,:);
    Be(4,ivdof) = dNi_dX(1,:);
    Be(5,ivdof) = dNi_dX(3,:);
    Be(5,iwdof) = dNi_dX(2,:);
    Be(6,iudof) = dNi_dX(3,:);
    Be(6,iwdof) = dNi_dX(1,:);
    
    SeGCSg(6*(iGauss-1)+(1:6),:) = C*Be*UeGCS; 
end

% for iGauss=1:nXi
%   Ni2 = sh_vs8(xi(iGauss,1),xi(iGauss,2),xi(iGauss,3));
%   Extrap2(iGauss,:) = Ni2; 
% end
% Extrap2 = kron(Extrap2,eye(6));

% extrapolate stress to nodes
xig = xi1D(1);
Extrap=[(xig - 1)*(xig - 1)*(xig - 1), (xig - 1)*(xig - 1)*(xig + 1), (xig - 1)*(xig + 1)*(xig - 1), (xig - 1)*(xig + 1)*(xig + 1), (xig + 1)*(xig - 1)*(xig - 1), (xig + 1)*(xig - 1)*(xig + 1), (xig + 1)*(xig + 1)*(xig - 1), (xig + 1)*(xig + 1)*(xig + 1)
(xig + 1)*(xig - 1)*(xig - 1), (xig + 1)*(xig - 1)*(xig + 1), (xig + 1)*(xig + 1)*(xig - 1), (xig + 1)*(xig + 1)*(xig + 1), (xig - 1)*(xig - 1)*(xig - 1), (xig - 1)*(xig - 1)*(xig + 1), (xig - 1)*(xig + 1)*(xig - 1), (xig - 1)*(xig + 1)*(xig + 1)
(xig + 1)*(xig + 1)*(xig - 1), (xig + 1)*(xig + 1)*(xig + 1), (xig + 1)*(xig - 1)*(xig - 1), (xig + 1)*(xig - 1)*(xig + 1), (xig - 1)*(xig + 1)*(xig - 1), (xig - 1)*(xig + 1)*(xig + 1), (xig - 1)*(xig - 1)*(xig - 1), (xig - 1)*(xig - 1)*(xig + 1)
(xig - 1)*(xig + 1)*(xig - 1), (xig - 1)*(xig + 1)*(xig + 1), (xig - 1)*(xig - 1)*(xig - 1), (xig - 1)*(xig - 1)*(xig + 1), (xig + 1)*(xig + 1)*(xig - 1), (xig + 1)*(xig + 1)*(xig + 1), (xig + 1)*(xig - 1)*(xig - 1), (xig + 1)*(xig - 1)*(xig + 1)
(xig - 1)*(xig - 1)*(xig + 1), (xig - 1)*(xig - 1)*(xig - 1), (xig - 1)*(xig + 1)*(xig + 1), (xig - 1)*(xig + 1)*(xig - 1), (xig + 1)*(xig - 1)*(xig + 1), (xig + 1)*(xig - 1)*(xig - 1), (xig + 1)*(xig + 1)*(xig + 1), (xig + 1)*(xig + 1)*(xig - 1)
(xig + 1)*(xig - 1)*(xig + 1), (xig + 1)*(xig - 1)*(xig - 1), (xig + 1)*(xig + 1)*(xig + 1), (xig + 1)*(xig + 1)*(xig - 1), (xig - 1)*(xig - 1)*(xig + 1), (xig - 1)*(xig - 1)*(xig - 1), (xig - 1)*(xig + 1)*(xig + 1), (xig - 1)*(xig + 1)*(xig - 1)
(xig + 1)*(xig + 1)*(xig + 1), (xig + 1)*(xig + 1)*(xig - 1), (xig + 1)*(xig - 1)*(xig + 1), (xig + 1)*(xig - 1)*(xig - 1), (xig - 1)*(xig + 1)*(xig + 1), (xig - 1)*(xig + 1)*(xig - 1), (xig - 1)*(xig - 1)*(xig + 1), (xig - 1)*(xig - 1)*(xig - 1)
(xig - 1)*(xig + 1)*(xig + 1), (xig - 1)*(xig + 1)*(xig - 1), (xig - 1)*(xig - 1)*(xig + 1), (xig - 1)*(xig - 1)*(xig - 1), (xig + 1)*(xig + 1)*(xig + 1), (xig + 1)*(xig + 1)*(xig - 1), (xig + 1)*(xig - 1)*(xig + 1), (xig + 1)*(xig - 1)*(xig - 1)];
Extrap = Extrap/(8*xig^3);
Extrap = kron(Extrap,eye(6));

SeGCS = Extrap*SeGCSg;
% SeGCS2 = Extrap2\SeGCSg;

% stresses in LCS
if nargout > 1
    t = trans_solid20(Node);
    theta = vtrans_solid(t);
    SeLCS = zeros(size(SeGCS)); 
    for iNode=1:8
    SeLCS((1:6)+6*(iNode-1),:)= theta*SeGCS((1:6)+6*(iNode-1),:);  
    end
    vLCS = t.'; vLCS = vLCS(:).';
end


if nargin>6 & ~isempty(gcs)
    
    switch lower(gcs)
        case {'cart'}
            
        case {'cyl'}
            a1 = [Node(1:8,1) Node(1:8,2) zeros(8,1)];
            a2 = [-Node(1:8,2) Node(1:8,1) zeros(8,1)];
            
            for iNode=1:8
                A= [a1(iNode,:)/norm(a1(iNode,:));a2(iNode,:)/norm(a2(iNode,:));0 0 1];
                theta = vtrans_solid(A);            
                SeGCS((1:6)+6*(iNode-1),:)= theta*SeGCS((1:6)+6*(iNode-1),:);
            end
            
        case {'sph'}
            a1 = Node(1:8,:);
            a2 = [-Node(1:8,2) Node(1:8,1) zeros(8,1)];
            a3 = [-Node(1:8,1).*Node(1:8,3) -Node(1:8,2).*Node(1:8,3) Node(1:8,1).^2+Node(1:8,2).^2];
            for iNode=1:8
                A= [a1(iNode,:)/norm(a1(iNode,:));a2(iNode,:)/norm(a2(iNode,:));a3(iNode,:)/norm(a3(iNode,:))];
                theta = vtrans_solid(A);
                SeGCS((1:6)+6*(iNode-1),:)= theta*SeGCS((1:6)+6*(iNode-1),:);
            end
            
    end
end


end

function z = cross(x,y)
z = x([2 1 1]).*y([3 3 2])-y([2 1 1]).*x([3 3 2]);
z(2) = -z(2);
end