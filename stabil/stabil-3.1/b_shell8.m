function [Bg,J] = b_shell8(Ni,dN_dxi,dN_deta,zeta,Node,h,v1i,v2i,v3i,d)

%B_SHELL8   b matrix for a shell8 element in global coordinate system.
%
%   [Bg,J] = b_shell8(Ni,dN_dxi,dN_deta,zetar,Node,h,v1i,v2i,v3i) returns
%   the element b matrix in the global coordinate system and the Jacobian
%   of the parametric transformation. Both are evaluated in the natural
%   coordinates (xi,eta and zetar) which were used to calculate
%   Ni,dN_dxi,dN_deta and zetar.
%
%   Node       Node definitions                             [x y z] (8 * 3)
%              Nodes should have the following order:
%              4----7----3
%              |         |
%              8         6
%              |         |
%              1----5----2
%   Ni         Shape functions for quadratic serendipity element    (8 * 1)
%   dN_dxi     first derivatives of shape functions Ni              (8 * 1)
%   dN_deta    first derivatives of shape functions Ni              (8 * 1)
%   h          scalar or vector containing thickness      scalar or (8 * 1)
%   v(1,2,3)i  components of the local coordinate system in node i  (8 * 3)
%   d          Nodal offset from shell mid plane       scalar (default = 0)
%   Bg         b matrix of shell8 element                          (6 * 48)
%   J          Jacobian of the parametric transformation            (3 * 3)
%
%   See also SE_SHELL8, KE_SHELL8.

% Miche Jansen
% 2009

% offset
if nargin==10
    zetan = 2*d./h;
else
    zetan = 0;
end
zetar = zeta-zetan;

a = Node(:,1)+zetar.*h.*v3i(:,1)/2;
b = Node(:,2)+zetar.*h.*v3i(:,2)/2;
c = Node(:,3)+zetar.*h.*v3i(:,3)/2;

%xg = Ni.'*a;
%yg = Ni.'*b;
%zg = Ni.'*c;

J = [dN_dxi.'*a dN_dxi.'*b dN_dxi.'*c;
    dN_deta.'*a dN_deta.'*b dN_deta.'*c;
    Ni.'*(h.*v3i(:,1)/2) Ni.'*(h.*v3i(:,2)/2) Ni.'*(h.*v3i(:,3)/2)];

invJ = inv(J);
% alle afgeleiden in het globale cartesische assenstelsel worden hier
% bepaald. De afgeleiden worden gesorteerd per DOF
du_dx_u = dN_dxi*invJ(1,1)+dN_deta*invJ(1,2);
%     dv_dx_v = du_dx_u;
%     dw_dx_w = du_dx_u;
%     du_dx_alfa = (du_dx_u).*zetar.*h.*v1i(:,1)/2+Ni.*h.*v1i(:,1)/2*invJ(1,3);
%     dv_dx_alfa = (du_dx_u).*zetar.*h.*v1i(:,2)/2+Ni.*h.*v1i(:,2)/2*invJ(1,3);
%     dw_dx_alfa = (du_dx_u).*zetar.*h.*v1i(:,3)/2+Ni.*h.*v1i(:,3)/2*invJ(1,3);
%     du_dx_beta = -(du_dx_u).*zetar.*h.*v2i(:,1)/2-Ni.*h.*v2i(:,1)/2*invJ(1,3);
%     dv_dx_beta = -(du_dx_u).*zetar.*h.*v2i(:,2)/2-Ni.*h.*v2i(:,2)/2*invJ(1,3);
%     dw_dx_beta = -(du_dx_u).*zetar.*h.*v2i(:,3)/2-Ni.*h.*v2i(:,3)/2*invJ(1,3);
%
du_dy_u = dN_dxi*invJ(2,1)+dN_deta*invJ(2,2);
%     dv_dy_v = du_dy_u;
%     dw_dy_w = du_dy_u;
%     du_dy_alfa = (du_dy_u).*zetar.*h.*v1i(:,1)/2+Ni.*h.*v1i(:,1)/2*invJ(2,3);
%     dv_dy_alfa = (du_dy_u).*zetar.*h.*v1i(:,2)/2+Ni.*h.*v1i(:,2)/2*invJ(2,3);
%     dw_dy_alfa = (du_dy_u).*zetar.*h.*v1i(:,3)/2+Ni.*h.*v1i(:,3)/2*invJ(2,3);
%     du_dy_beta = -(du_dy_u).*zetar.*h.*v2i(:,1)/2-Ni.*h.*v2i(:,1)/2*invJ(2,3);
%     dv_dy_beta = -(du_dy_u).*zetar.*h.*v2i(:,2)/2-Ni.*h.*v2i(:,2)/2*invJ(2,3);
%     dw_dy_beta = -(du_dy_u).*zetar.*h.*v2i(:,3)/2-Ni.*h.*v2i(:,3)/2*invJ(2,3);
%
du_dz_u = dN_dxi*invJ(3,1)+dN_deta*invJ(3,2);
%     dv_dz_v = du_dz_u;
%     dw_dz_w = du_dz_u;
%     du_dz_alfa = (du_dz_u).*zetar.*h.*v1i(:,1)/2+Ni.*h.*v1i(:,1)/2*invJ(3,3);
%     dv_dz_alfa = (du_dz_u).*zetar.*h.*v1i(:,2)/2+Ni.*h.*v1i(:,2)/2*invJ(3,3);
%     dw_dz_alfa = (du_dz_u).*zetar.*h.*v1i(:,3)/2+Ni.*h.*v1i(:,3)/2*invJ(3,3);
%     du_dz_beta = -(du_dz_u).*zetar.*h.*v2i(:,1)/2-Ni.*h.*v2i(:,1)/2*invJ(3,3);
%     dv_dz_beta = -(du_dz_u).*zetar.*h.*v2i(:,2)/2-Ni.*h.*v2i(:,2)/2*invJ(3,3);
%     dw_dz_beta = -(du_dz_u).*zetar.*h.*v2i(:,3)/2-Ni.*h.*v2i(:,3)/2*invJ(3,3);


Bg = zeros(6,48);
Bg(1,1:6:43) = du_dx_u;
Bg(1,4:6:46) = du_dx_u.*zetar.*h.*v1i(:,1)/2+Ni.*h.*v1i(:,1)/2*invJ(1,3);
Bg(1,5:6:47) = -du_dx_u.*zetar.*h.*v2i(:,1)/2-Ni.*h.*v2i(:,1)/2*invJ(1,3);
Bg(2,2:6:44) = du_dy_u;
Bg(2,4:6:46) = du_dy_u.*zetar.*h.*v1i(:,2)/2+Ni.*h.*v1i(:,2)/2*invJ(2,3);
Bg(2,5:6:47) = -du_dy_u.*zetar.*h.*v2i(:,2)/2-Ni.*h.*v2i(:,2)/2*invJ(2,3);
Bg(3,3:6:45) = du_dz_u;
Bg(3,4:6:46) = du_dz_u.*zetar.*h.*v1i(:,3)/2+Ni.*h.*v1i(:,3)/2*invJ(3,3);
Bg(3,5:6:47) = -du_dz_u.*zetar.*h.*v2i(:,3)/2-Ni.*h.*v2i(:,3)/2*invJ(3,3);
Bg(4,1:6:43) = du_dy_u;
Bg(4,2:6:44) = du_dx_u;
Bg(4,4:6:46) = du_dy_u.*zetar.*h.*v1i(:,1)/2+Ni.*h.*v1i(:,1)/2*invJ(2,3) + du_dx_u.*zetar.*h.*v1i(:,2)/2+Ni.*h.*v1i(:,2)/2*invJ(1,3);
Bg(4,5:6:47) = -du_dy_u.*zetar.*h.*v2i(:,1)/2-Ni.*h.*v2i(:,1)/2*invJ(2,3)-du_dx_u.*zetar.*h.*v2i(:,2)/2-Ni.*h.*v2i(:,2)/2*invJ(1,3);
Bg(5,2:6:44) = du_dz_u;
Bg(5,3:6:45) = du_dy_u;
Bg(5,4:6:46) = du_dz_u.*zetar.*h.*v1i(:,2)/2+Ni.*h.*v1i(:,2)/2*invJ(3,3)+du_dy_u.*zetar.*h.*v1i(:,3)/2+Ni.*h.*v1i(:,3)/2*invJ(2,3);
Bg(5,5:6:47) = -du_dz_u.*zetar.*h.*v2i(:,2)/2-Ni.*h.*v2i(:,2)/2*invJ(3,3)-du_dy_u.*zetar.*h.*v2i(:,3)/2-Ni.*h.*v2i(:,3)/2*invJ(2,3);
Bg(6,1:6:43) = du_dz_u;
Bg(6,3:6:45) = du_dx_u;
Bg(6,4:6:46) = du_dz_u.*zetar.*h.*v1i(:,1)/2+Ni.*h.*v1i(:,1)/2*invJ(3,3)+du_dx_u.*zetar.*h.*v1i(:,3)/2+Ni.*h.*v1i(:,3)/2*invJ(1,3);
Bg(6,5:6:47) = -du_dz_u.*zetar.*h.*v2i(:,1)/2-Ni.*h.*v2i(:,1)/2*invJ(3,3)-du_dx_u.*zetar.*h.*v2i(:,3)/2-Ni.*h.*v2i(:,3)/2*invJ(1,3);


end

