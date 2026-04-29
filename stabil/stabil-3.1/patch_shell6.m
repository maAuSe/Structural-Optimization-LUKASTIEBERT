function [pxyz,pind,pvalue]=patch_shell6(Nodes,NodeNum,Values)

%PATCH_SHELL6  Patch information of the shell6 elements for plotting.
%
%   [pxyz,pind,pvalue] = patch_shell6(Nodes,NodeNum,Values) returns matrices
%   to plot patches of shell6 elements.
%
%   Nodes      Node definitions        [NodID x y z]
%   NodeNum    Node numbers       [NodID1 NodID2 NodID3 NodID4] (nElem * 6)
%   Values     Values assigned to nodes used for coloring    (nElem * 6)
%   pxyz       Coordinates of Nodes                          (6*nElem * 3)
%   pind       indices of Nodes                                (nElem * 8)
%   pvalue     Values arranged per Node                      (6*nElem * 1)
%
%   See also PLOTSTRESSCONTOURF, PLOTSHELLFCONTOURF.

% Miche Jansen
% 2010

nElem=size(NodeNum,1);
nNode = 6;
NodeNum = NodeNum(:,1:nNode).';
nPoints = 5;

% check nodes
if length(unique(Nodes(:,1)))~=length(Nodes(:,1))
    for iNode = 1:length(Nodes(:,1))
        loc = find(abs(Nodes(iNode,1)-Nodes(:,1))<1e-4);
        if length(loc)>1
        error('Node %i is multiply defined.',Nodes(iNode,1))
        end
    end
end

Xn=Nodes(:,2);
Yn=Nodes(:,3);
Zn=Nodes(:,4);


N = zeros(4*(nPoints),6);
s = linspace(0,1,nPoints).';
st = [s zeros(nPoints,1);
      1-s s;
      zeros(nPoints,1) flipud(s)];

for ind = 1:3*(nPoints)
    Ni = sh_t6(st(ind,1),st(ind,2));
    N(ind,:) = Ni; 
end

if size(Values,2) < 3
Values = Values(:,ones(1,6));    
elseif size(Values,2) < 6 
    Values = [Values(:,1),Values(:,2),Values(:,3),(Values(:,1)+Values(:,2))/2,...
              (Values(:,2)+Values(:,3))/2,(Values(:,3)+Values(:,1))/2];
else 
Values = Values(:,1:6); 
end

Values=Values.';

[exis,loc] = ismember(round(NodeNum),round(Nodes(:,1)));

if ~all(exis(:))
   unknowns = NodeNum(~exis);
   error('Node %i is not defined.',unknowns(1));
end

Xn = Xn(loc);
Yn = Yn(loc);
Zn = Zn(loc);

px=N*Xn;
py=N*Yn;
pz=N*Zn;

pvalue=N*Values;

pind = reshape((1:size(N,1)*nElem),size(N,1),nElem);
pind = pind.';

pxyz = [px(:),py(:),pz(:)];

pvalue = pvalue(:);
end