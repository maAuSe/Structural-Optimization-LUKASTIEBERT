function [BeGCS] = be_plane3(Node,Section,Material,UeGCS,Options,gcs)

x=Node(:,1);
y=Node(:,2);

A12=1/2*(x(1)*y(2)-x(2)*y(1));
A23=1/2*(x(2)*y(3)-x(3)*y(2));
A31=1/2*(x(3)*y(1)-x(1)*y(3));

A=A12+A23+A31;

x32=x(3)-x(2);
x13=x(1)-x(3);
x21=x(2)-x(1);

y23=y(2)-y(3);
y31=y(3)-y(1);
y12=y(1)-y(2);

B=1/2/A*[y23   0  y31   0  y12   0
         0   x32    0 x13    0 x21
         x32 y23  x13 y31  x21 y12];

strain=B*UeGCS;     
         
BeGCS=repmat(strain,3,1);