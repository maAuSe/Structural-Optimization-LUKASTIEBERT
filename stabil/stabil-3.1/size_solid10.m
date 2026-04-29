function [S] = size_solid10(Node)
% size_solid10   Compute solid10 element size (volume).            
%                                                                            
%     S = size_solid10(NodeNum) computes the size of a solid10 element.                                                       
%                                                                            
%     Node        Node definitions                    [x y z] (10 * 3)
%     S           Element size


% ! Constant metric !

x=Node(:,1);
y=Node(:,2);
z=Node(:,3);

V01=1/6*(x(2)*y(3)*z(4)-x(2)*y(4)*z(3)-x(3)*y(2)*z(4)+x(3)*y(4)*z(2)+x(4)*y(2)*z(3)-x(4)*y(3)*z(2));
V02=1/6*(x(1)*y(4)*z(3)-x(1)*y(3)*z(4)+x(3)*y(1)*z(4)-x(3)*y(4)*z(1)-x(4)*y(1)*z(3)+x(4)*y(3)*z(1));
V03=1/6*(x(1)*y(2)*z(4)-x(1)*y(4)*z(2)-x(2)*y(1)*z(4)+x(2)*y(4)*z(1)+x(4)*y(1)*z(2)-x(4)*y(2)*z(1));
V04=1/6*(x(1)*y(3)*z(2)-x(1)*y(2)*z(3)+x(2)*y(1)*z(3)-x(2)*y(3)*z(1)-x(3)*y(1)*z(2)+x(3)*y(2)*z(1));

S=V01+V02+V03+V04;  % element volume
