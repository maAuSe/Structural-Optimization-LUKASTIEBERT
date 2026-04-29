function [s,t,NodeNum,Elements] = grid_shell4(m,n,Type,Section,Material)

%GRID_SHELL4  Grid in natural coordinates for mapped meshing.
%
%   [s,t,NodeNum,Elements] = grid_shell4(m,n,Type,Section,Material)
%   returns matrices of a grid in the natural coordinate system, which can 
%   be used for mapped meshing.
%
%   s           s-coordinate of nodes  (1 * nNodes)
%   t           t-coordinate of nodes  (1 * nNodes)
%   NodeNum     Node numbers order on grid ((m+1) * (n+1))        
%   Elements    Node numbers are saved per element here (nElem * 8)
%   Type        Type ID of meshed elements   
%   Section     Section ID of meshed elements
%   Material    Material ID of meshed elements
%
%   See also MAKEMESH, GRID_SHELL8.

% Miche Jansen
% 2009

s = linspace(-1,1,n+1);
t = fliplr(linspace(-1,1,m+1));

nNodes = (m+1)*(n+1);
NodeNum = reshape((1:nNodes),m+1,n+1);

nElem = m*n;
Elements = zeros(nElem,8);

for j = 1:n
    for i = 1:m
        Elements(i+(j-1)*m,:) = [i+(j-1)*m cell2mat(Type(1)) Section(1) Material(1) NodeNum(i+1,j) NodeNum(i+1,j+1) NodeNum(i,j+1) NodeNum(i,j)];
    end
end
end