function [F,dFdx] = loads_truss(DLoad,Node,dDLoaddx,dNodedx)

%LOADS_TRUSS   Equivalent nodal forces for a truss element in the GCS.
%
%   F = LOADS_TRUSS(DLoad,Node)
%   computes the equivalent nodal forces of a distributed load 
%   (in the global coordinate system).
%
%   [F,dFdx] = LOADS_TRUSS(DLoad,Node,dDLoaddx,dNodedx)
%   additionally computes the derivatives of the equivalent nodal forces
%   with respect to the design variables x.
%
%   DLoad      Distributed load          [n1globalX n1globalY n1globalZ ...]
%                   (6/8 * nLC)
%   Node       Node definitions          [x y z] (2 * 3)
%   dDLoaddx   Distributed load derivatives     (6/8 * nLC * nVar)
%   dNodedx    Node definitions derivatives     (SIZE(Node) * nVar)
%   F          Load vector              (6 * nLC)
%   dFdx       Load vector derivatives 	(6 * nLC * nVar)
%
%   See also ELEMLOADS, LOADS_BEAM.

% David Dooms, Wouter Dillen
% October 2008, April 2017

if nargin<3, dDLoaddx = []; end
if nargin<4, dNodedx = []; end

nVar = 0;
if nargout>1 && (~isempty(dDLoaddx) || ~isempty(dNodedx))
    nVar = max([size(dDLoaddx,3),size(dNodedx,3)]);
end 

if nVar==0 || isempty(dDLoaddx), dDLoaddx = zeros(size(DLoad,1),size(DLoad,2),nVar,size(DLoad,3)); end
if nVar==0 || isempty(dNodedx), dNodedx = zeros([size(Node),nVar]); end


% check input
if size(DLoad,3)>1
    error('Multiple distributed loads on a single truss element.');
elseif size(DLoad,1)==8 
    if any(any(DLoad(7:8,:)))~=0
        error('Partial distributed loads are not allowed on truss elements.');
    end
end

% Element length
L=norm(Node(2,:)-Node(1,:));
dLdx=(permute(dNodedx(2,:,:)-dNodedx(1,:,:),[3 2 1])*(Node(2,:)-Node(1,:)).')/L;

% F
F=zeros(6,size(DLoad,2));
F(1:3,:)=L/6*(2*DLoad(1:3,:)+DLoad(4:6,:));
F(4:6,:)=L/6*(DLoad(1:3,:)+2*DLoad(4:6,:));
dFdx=zeros(size(F,1),size(F,2),nVar);
for n=1:nVar
    dFdx(1:3,:,n)=F(1:3,:)/L*dLdx(n) + L/6*(2*dDLoaddx(1:3,:,n)+dDLoaddx(4:6,:,n));
    dFdx(4:6,:,n)=F(4:6,:)/L*dLdx(n) + L/6*(dDLoaddx(1:3,:,n)+2*dDLoaddx(4:6,:,n));
end
