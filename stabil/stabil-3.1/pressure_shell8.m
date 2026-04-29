function F = pressure_shell8(Pressure,Node)

%PRESSURE_SHELL8   Equivalent nodal forces for a shell8 element in the GCS
%                  due to a pressure normal to the element surface.
%
%   F = pressure_shell(Pressure,Node)
%   computes the equivalent nodal forces of a pressure load normal to 
%   the elements surface.   
%
%   Pressure   Distributed loads in corner Nodes    [p1localZ; ...]
%                                                                   (4 * 1)
%   Node       Node definitions           [x y z] (8 * 3)
%   F          Load vector  (48 * 1)
%
%   See also  ELEMPRESSURE, PRESSURE_SHELL4.

% Miche Jansen
% 2009

Node = Node(1:8,1:3);
% integratiepunten
[x,H] = gaussq(2);

F = zeros(48,size(Pressure,2));


xi_eta = [ -1 -1;1  -1;1 1; -1 1;0 -1; 1 0; 0 1; -1 0];
v3i = zeros(8,3);
for iNode = 1:8
    [Ni,dN_dxi,dN_deta] = sh_qs8(xi_eta(iNode,1),xi_eta(iNode,2));
    Jm = [dN_dxi.'*Node(:,1) dN_dxi.'*Node(:,2) dN_dxi.'*Node(:,3);
         dN_deta.'*Node(:,1) dN_deta.'*Node(:,2) dN_deta.'*Node(:,3)];
    v3i(iNode,:) = cross(Jm(1,:),Jm(2,:));
    v3i(iNode,:) = v3i(iNode,:)/norm(v3i(iNode,:));
end  

Pressure = [Pressure;(Pressure(1,:)+Pressure(2,:))/2;
            (Pressure(2,:)+Pressure(3,:))/2;
            (Pressure(3,:)+Pressure(4,:))/2;
            (Pressure(4,:)+Pressure(1,:))/2];
%p148 Zienkiewicz deel1 (2005)
for iGauss=1:size(x,1)
    xi = x(iGauss,1);
    eta = x(iGauss,2);
    [Ni,dN_dxi,dN_deta] = sh_qs8(xi,eta);
    N = zeros(3,48);
    N(1,1:6:43) = Ni;
    N(2,2:6:44) = Ni;
    N(3,3:6:45) = Ni;
    Pressureg = Ni.'*Pressure;
    v3g = Ni.'*v3i;
    v3g = v3g.'/norm(v3g);
    det = norm(cross(Node.'*dN_dxi,Node.'*dN_deta));

    F = F + H(iGauss)*N.'*v3g*Pressureg*det;
end