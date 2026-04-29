function F = pressure_shell4(Pressure,Node)

%PRESSURE_SHELL4   Equivalent nodal forces for a shell4 element in the GCS
%                  due to a pressure normal to the element surface.
%
%   F = pressure_shell4(Pressure,Node)
%   computes the equivalent nodal forces of a pressure load normal to the 
%   elements surface. 
%
%   Pressure   Distributed pressure in corner Nodes    [p1lobalZ; ...]
%                                                                   (4 * 1)
%   Node       Node definitions           [x y z] (4 * 3)
%   F          Load vector  (24 * 1)
%
%   See also  ELEMPRESSURE, PRESSURE_SHELL8.

% Miche Jansen
% 2009

Node = Node(1:4,1:3);
[t,Node_lc]=trans_shell4(Node);

FLCS = zeros(24,size(Pressure,2));

seq=[1 2 3;
     1 3 4;
     1 2 4;
     2 3 4];
 
x = 1/6*[4 1 1;1 4 1;1 1 4];
H = [1 1 1]/3;
for n =1:4
    for iGauss=1:3
        Node_t=Node_lc(seq(n,:),:);
        Pressureg = x(iGauss,:)*Pressure(seq(n,:),:);
        ind = 6*([seq(n,1)*ones(3,1);seq(n,2)*ones(3,1);seq(n,3)*ones(3,1)]-1)+repmat([3; 4; 5],3,1);
        b = [Node_t(2,2)-Node_t(3,2);
             Node_t(3,2)-Node_t(1,2);
             Node_t(1,2)-Node_t(2,2)];
        c = [Node_t(3,1)-Node_t(2,1);
             Node_t(1,1)-Node_t(3,1);
             Node_t(2,1)-Node_t(1,1)];
        Det = b(2)*c(3)-b(3)*c(2);
        Ni = sh_t(x(:,iGauss),b,c);
        FLCS(ind,:) = FLCS(ind,:) +H(iGauss)*Ni.'*Pressureg*Det/4;
    end
end


T=blkdiag(t,t,t,t,t,t,t,t);

F=T.'*FLCS; 

end