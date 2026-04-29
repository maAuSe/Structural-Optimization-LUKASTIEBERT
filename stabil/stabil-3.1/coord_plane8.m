function [X,Y,Z]=coord_plane8(Nodes,NodeNum)

%COORD_PLANE8  Coordinates of the shell8 elements for plotting.
%
%   coord_plane8(Nodes,NodeNum)
%   returns the coordinates of the plane8 elements for plotting.
%
%   Nodes      Node definitions        [NodID x y z] (nNodes * 4)
%   NodeNum    Node numbers            [NodID1 NodID2 NodID3 NodID4] (nElem * 8)
%   X          X coordinates  (20 * nElem)    
%   Y          Y coordinates  (20 * nElem)
%   Z          Z coordinates  (20 * nElem)
%
%   See also COORD_TRUSS, PLOTELEM.

% Miche Jansen
% 2009

if ~isempty(NodeNum)
nElem=size(NodeNum,1);
nNode = size(NodeNum,2);
nPoints = 5; % aanpassen voor betere weergave (indien element gekromd is)(formules onderaan uncommenten indien enkel knopen juist moeten zijn)

Xn=zeros(nNode,nElem);
Yn=zeros(nNode,nElem);
Zn=zeros(nNode,nElem);

N = zeros(4*(nPoints),8);
s = linspace(-1,1,nPoints).';
st = [s -ones(nPoints,1);
      ones(nPoints,1) s;
      -s ones(nPoints,1);
      -ones(nPoints,1) -s];


for ind = 1:4*(nPoints)
    Ni = sh_qs8(st(ind,1),st(ind,2));
    N(ind,:) = Ni; 
end

for iElem=1:nElem
    for ind =1:8 
    loc=find(Nodes(:,1)==NodeNum(iElem,ind));
    if isempty(loc)
        error('Node %i is not defined.',NodeNum(iElem,1))
    elseif length(loc)>1
        error('Node %i is multiply defined.',NodeNum(iElem,1))
    else
        Xn(ind,iElem)=Nodes(loc,2);
        Yn(ind,iElem)=Nodes(loc,3);
        Zn(ind,iElem)=Nodes(loc,4);
    end
    end
   
end
X = N*Xn;
Y = N*Yn;
Z = N*Zn;

else X = [];Y = []; Z = [];
end

% X = [Xn(1,:);Xn(5,:);Xn(2,:);Xn(6,:);Xn(3,:);Xn(7,:);Xn(4,:);Xn(8,:)];
% Y = [Yn(1,:);Yn(5,:);Yn(2,:);Yn(6,:);Yn(3,:);Yn(7,:);Yn(4,:);Yn(8,:)];
% Z = [Zn(1,:);Zn(5,:);Zn(2,:);Zn(6,:);Zn(3,:);Zn(7,:);Zn(4,:);Zn(8,:)];

end