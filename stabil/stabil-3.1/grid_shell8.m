function [s,t,NodeNum,Elements] = grid_shell8(m,n,Type,Section,Material)

%GRID_SHELL8  Grid in natural coordinates for mapped meshing.
%
%   [s,t,NodeNum,Elements] = grid_shell8(m,n,Type,Section,Material)
%   returns matrices of a grid in the natural coordinate system, which can 
%   be used for mapped meshing.
%
%   s           s-coordinate of nodes  (1 * nNodes)
%   t           t-coordinate of nodes  (1 * nNodes)
%   NodeNum     Node numbers order on grid ((m+1) * (n+1))        
%   Elements    Node numbers are saved per element here (nElem * 12)
%   Type        Type ID of meshed elements   
%   Section     Section ID of meshed elements
%   Material    Material ID of meshed elements
%
%   See also MAKEMESH, GRID_SHELL4.

% Miche Jansen
% 2009

s = linspace(-1,1,2*n+1);
t = fliplr(linspace(-1,1,2*m+1));

NodeNum = zeros(2*m+1,2*n+1);
for j=1:n
    for i=1:m
        NodeNum([2*(i-1)+1:2*(i-1)+3],[2*(j-1)+1:2*(j-1)+3]) = [2*(i-1)+1+(3*m+2)*(j-1) i+(2*m+1)+(3*m+2)*(j-1) 2*(i-1)+1+(3*m+2)*j;
                                        2*(i-1)+2+(3*m+2)*(j-1) NaN 2*(i-1)+2+(3*m+2)*j;
                                        2*(i-1)+3+(3*m+2)*(j-1) i+1+(2*m+1)+(3*m+2)*(j-1) 2*(i-1)+3+(3*m+2)*j];
    end
end


nElem = m*n;
Elements = zeros(nElem,12);

for j = 1:n
    for i = 1:m
        Elements(i+(j-1)*m,:) = [i+(j-1)*(m) cell2mat(Type(1)) Section(1) Material(1) NodeNum(2*(i-1)+3,2*(j-1)+1) NodeNum(2*(i-1)+3,2*(j-1)+3) NodeNum(2*(i-1)+1,2*(j-1)+3) NodeNum(2*(i-1)+1,2*(j-1)+1), ...
             NodeNum(2*(i-1)+3,2*(j-1)+2) NodeNum(2*(i-1)+2,2*(j-1)+3) NodeNum(2*(i-1)+1,2*(j-1)+2) NodeNum(2*(i-1)+2,2*(j-1)+1)];
    end
end

end
