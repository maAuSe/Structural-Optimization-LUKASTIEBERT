function [S] = size_plane6(Node)
% size_plane6   Compute plane6 element size (area).              
%                                                                            
%     s = size_plane6(NodeNum) computes the size (area) of a plane6 element.                                                       
%                                                                            
%     Node        Node definitions                    [x y z] (6 * 3)
%     S           Element size
%


% ! Constant metric !

x=Node(:,1);
y=Node(:,2);

A12=1/2*(x(1)*y(2)-x(2)*y(1));
A23=1/2*(x(2)*y(3)-x(3)*y(2));
A31=1/2*(x(3)*y(1)-x(1)*y(3));

S=A12+A23+A31; % Element area