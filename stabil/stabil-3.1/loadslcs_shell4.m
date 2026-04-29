function FLCS = loadslcs_shell4(DLoadLCS,Node_lc)

%LOADSLCS_SHELL4   Equivalent nodal forces for a shell4 element in the LCS.
%
%   F = loadslcs_shell4(DLoadLCS,Node)
%   computes the equivalent nodal forces of a distributed load 
%   (in the local coordinate system).
%
%   DLoadLCS   Distributed loads      [n1localX; n1localY; n1localZ; ...]
%              in corner Nodes         (12 * 1)
%   Node       Node definitions       [x y z] (4 * 3)
%   FLCS       Load vector  (24 * 1)
%
%   See also LOADSLCS_BEAM, ELEMLOADS, LOADS_TRUSS.

% Miche Jansen
% 2009

% integratiepunten
[x,H] = gaussq(2);

FLCS = zeros(24,size(DLoadLCS,2));

for iGauss=1:4
    xi = x(iGauss,1);
    eta = x(iGauss,2);
    [Ni,dN_dxi,dN_deta] = sh_qs4(xi,eta);
    N = zeros(2,24);
    N(1,1:6:19) = Ni;
    N(2,2:6:20) = Ni;
    
    J = [dN_dxi.'*Node_lc(:,1) dN_dxi.'*Node_lc(:,2);
         dN_deta.'*Node_lc(:,1) dN_deta.'*Node_lc(:,2)];
     
    Dloadg = [Ni.'*DLoadLCS(1:3:10,:); Ni.'*DLoadLCS(2:3:11,:)];
    FLCS = FLCS + H(iGauss)*N.'*Dloadg*det(J);
end

seq=[1 2 3;
     1 3 4;
     1 2 4;
     2 3 4];
 
%x = 1/6*[4 1 1;1 4 1;1 1 4];
x = [1 0 0;0 1 0;0 0 1];
H = [1 1 1]/3;
for n =1:4
    for iGauss=1:3
        Node_t=Node_lc(seq(n,:),:);
        Dloadg = x(iGauss,:)*DLoadLCS(3*seq(n,:),:);
        ind = 6*([seq(n,1)*ones(3,1);seq(n,2)*ones(3,1);seq(n,3)*ones(3,1)]-1)+repmat([3; 4; 5],3,1);
        b = [Node_t(2,2)-Node_t(3,2);
             Node_t(3,2)-Node_t(1,2);
             Node_t(1,2)-Node_t(2,2)];
        c = [Node_t(3,1)-Node_t(2,1);
             Node_t(1,1)-Node_t(3,1);
             Node_t(2,1)-Node_t(1,1)];
        Det = b(2)*c(3)-b(3)*c(2);
        Ni = sh_t(x(:,iGauss),b,c);
        FLCS(ind,:) = FLCS(ind,:) +H(iGauss)*Ni.'*Dloadg*Det/4;
    end
end

end