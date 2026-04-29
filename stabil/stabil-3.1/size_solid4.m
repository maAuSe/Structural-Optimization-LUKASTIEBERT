function [S] = size_solid4(Node)
% size_solid4   Compute solid4 element size (volume).
%                                                                            
%     s = size_solid4(NodeNum) computes the size (volume) of a solid4 element.                                                       
%                                                                            
%     Node        Node definitions                    [x y z] (4 * 3)
%     S           Element size

x=Node(:,1);
y=Node(:,2);
z=Node(:,3);

V01=1/6*(x(2)*y(3)*z(4)-x(2)*y(4)*z(3)-x(3)*y(2)*z(4)+x(3)*y(4)*z(2)+x(4)*y(2)*z(3)-x(4)*y(3)*z(2));
V02=1/6*(x(1)*y(4)*z(3)-x(1)*y(3)*z(4)+x(3)*y(1)*z(4)-x(3)*y(4)*z(1)-x(4)*y(1)*z(3)+x(4)*y(3)*z(1));
V03=1/6*(x(1)*y(2)*z(4)-x(1)*y(4)*z(2)-x(2)*y(1)*z(4)+x(2)*y(4)*z(1)+x(4)*y(1)*z(2)-x(4)*y(2)*z(1));
V04=1/6*(x(1)*y(3)*z(2)-x(1)*y(2)*z(3)+x(2)*y(1)*z(3)-x(2)*y(3)*z(1)-x(3)*y(1)*z(2)+x(3)*y(2)*z(1));

s=V01+V02+V03+V04;  % element volume
