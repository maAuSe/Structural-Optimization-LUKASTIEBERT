function [t,Node_lc,W] = trans_shell4(Node,varargin)

%TRANS_SHELL4   Transform coordinate system for a shell4 element.
%
%   [t,Node_lc,W] = trans_shell4(Node)
%   [t,Node_lc]   = trans_shell4(Node)
%    t            = trans_shell4(Node)
%   computes the transformation matrix between the local and the global
%   coordinate system and the correction matrix for non-coplanar nodes
%   for the shell4 element.
%
%   Node       Node definitions           [x y z] (4 * 3)
%   t          Transformation matrix  (3 * 3)
%   Node_lc    Nodes in LCS               [x y z] (4 * 3)
%   W          Correction matrix for warped elements
%                                     (24 * 24)  
%
%   See also KE_BEAM, TRANS_TRUSS.

% Miche Jansen
% 2009

Nx=Node(2,:)-Node(1,:);        
Nx=Nx/norm(Nx);                

V24=Node(4,:)-Node(2,:);       % vector between nodes 2 and 4
V13=Node(3,:)-Node(1,:);       % vector between nodes 1 and 3

Nz=cross(V13,V24);             % vector along the local z-axis
if Nz==0
    error('Three nodes of the shell element are collinear')
end
Nz=Nz/norm(Nz);                % normalized vector along the local z-axis
Ny=cross(Nz,Nx);               
Nx=cross(Ny,Nz);               
t=[Nx/norm(Nx); Ny/norm(Ny); Nz];


ecenter = mean(Node,1);
Node_lc = (Node-ecenter([1 1 1 1],:))*t.';



% Warpingf = 2*abs(Node_lc(1,3))/h;
% if Warpingf > 0.1
%     warning('The shell4 element is strongly warped which could lead to inaccurate solutions (warping factor=%d)',Warpingf)
% end
W = eye(24);
if any(Node_lc(:,3) ~=0)
n1 = cross(Node_lc(2,:)-Node_lc(1,:),Node_lc(4,:)-Node_lc(1,:));
n2 = cross(Node_lc(3,:)-Node_lc(2,:),Node_lc(1,:)-Node_lc(2,:));   
n3 = cross(Node_lc(4,:)-Node_lc(3,:),Node_lc(2,:)-Node_lc(3,:));     
n4 = cross(Node_lc(1,:)-Node_lc(4,:),Node_lc(3,:)-Node_lc(4,:));         

n1 = n1/norm(n1);
n2 = n2/norm(n2);
n3 = n3/norm(n3);
n4 = n4/norm(n4);

fac = (Node_lc(1,1)*Node_lc(2,2)-Node_lc(2,1)*Node_lc(1,2)-Node_lc(1,1)*Node_lc(4,2)...
       + Node_lc(2,1)*Node_lc(3,2)-Node_lc(3,1)*Node_lc(2,2)+Node_lc(4,1)*Node_lc(1,2)...
       + Node_lc(3,1)*Node_lc(4,2)-Node_lc(4,1)*Node_lc(3,2));


% force correction
W(3,(1:6:19))=Node_lc(1,3)*(Node_lc(2,2)-Node_lc(4,2))*[1 -1 1 -1]/fac;
W(3,(2:6:20))=-Node_lc(1,3)*(Node_lc(2,1)-Node_lc(4,1))*[1 -1 1 -1]/fac;
W(9,(1:6:19))=-Node_lc(1,3)*(Node_lc(1,2)-Node_lc(3,2))*[1 -1 1 -1]/fac;
W(9,(2:6:20))=Node_lc(1,3)*(Node_lc(1,1)-Node_lc(3,1))*[1 -1 1 -1]/fac;    
W(15,(1:6:19))=-Node_lc(1,3)*(Node_lc(2,2)-Node_lc(4,2))*[1 -1 1 -1]/fac;
W(15,(2:6:20))=Node_lc(1,3)*(Node_lc(2,1)-Node_lc(4,1))*[1 -1 1 -1]/fac;
W(21,(1:6:19))=Node_lc(1,3)*(Node_lc(1,2)-Node_lc(3,2))*[1 -1 1 -1]/fac;
W(21,(2:6:20))=-Node_lc(1,3)*(Node_lc(1,1)-Node_lc(3,1))*[1 -1 1 -1]/fac;

% de volgende regels zijn enkel nodig indien de virtuele arbeid op de
% volledige Fz worden toegepast
facz1= -((Node_lc(2,1)*Node_lc(3,2))-(Node_lc(3,1)*Node_lc(2,2))-(Node_lc(2,1)*Node_lc(4,2))+(Node_lc(4,1)*Node_lc(2,2))+(Node_lc(3,1)*Node_lc(4,2))-(Node_lc(4,1)*Node_lc(3,2)))/2;
facz2= +((Node_lc(1,1)*Node_lc(3,2))-(Node_lc(3,1)*Node_lc(1,2))-(Node_lc(1,1)*Node_lc(4,2))+(Node_lc(4,1)*Node_lc(1,2))+(Node_lc(3,1)*Node_lc(4,2))-(Node_lc(4,1)*Node_lc(3,2)))/2;
facz3= -((Node_lc(1,1)*Node_lc(2,2))-(Node_lc(2,1)*Node_lc(1,2))-(Node_lc(1,1)*Node_lc(4,2))+(Node_lc(4,1)*Node_lc(1,2))+(Node_lc(2,1)*Node_lc(4,2))-(Node_lc(4,1)*Node_lc(2,2)))/2;
facz4= +((Node_lc(1,1)*Node_lc(2,2))-(Node_lc(2,1)*Node_lc(1,2))-(Node_lc(1,1)*Node_lc(3,2))+(Node_lc(3,1)*Node_lc(1,2))+(Node_lc(2,1)*Node_lc(3,2))-(Node_lc(3,1)*Node_lc(2,2)))/2;

W(3,(3:6:21))=W(3,(3:6:21))+Node_lc(1,3)*[1 -1 1 -1]*facz1/fac;
W(9,(3:6:21))=W(9,(3:6:21))+Node_lc(1,3)*[1 -1 1 -1]*facz2/fac;
W(15,(3:6:21))=W(15,(3:6:21))+Node_lc(1,3)*[1 -1 1 -1]*facz3/fac;
W(21,(3:6:21))=W(21,(3:6:21))+Node_lc(1,3)*[1 -1 1 -1]*facz4/fac;


%moment correction
W(6,4)=-n1(1)/n1(3); 
W(6,5)=-n1(2)/n1(3);
W(12,10)=-n2(1)/n2(3);
W(12,11)=-n2(2)/n2(3);
W(18,16)=-n3(1)/n3(3);
W(18,17)=-n3(2)/n3(3);
W(24,22)=-n4(1)/n4(3);
W(24,23)=-n4(2)/n4(3);

fac2 = -Node_lc(1,1)-Node_lc(1,2)+Node_lc(2,1)-Node_lc(2,2)...
       +Node_lc(3,1)+Node_lc(3,2)-Node_lc(4,1)+Node_lc(4,2);
   
R = -[W(6,4) W(6,5) W(12,10) W(12,11) W(18,16) W(18,17) W(24,22) W(24,23)]/fac2;

W([1 2 7 8 13 14 19 20],[4 5 10 11 16 17 22 23])=[R;-R;R;R;-R;R;-R;-R];
      
end
end

function z = cross(x,y)
z = x([2 1 1]).*y([3 3 2])-y([2 1 1]).*x([3 3 2]);
z(2) = -z(2);
end