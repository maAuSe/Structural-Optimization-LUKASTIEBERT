function [t,dtdx] = trans_beam(Node,dNodedx)

%TRANS_BEAM   Transform coordinate system for a beam element.
%
%   t = TRANS_BEAM(Node)
%       computes the transformation matrix between the local and the global
%       coordinate system for the beam element.
%   [t,dtdx] = TRANS_BEAM(Node,dNodedx)
%       additionally computes the derivatives of the transformation matrix
%       with respect to the design variables x.
%
%   Node       Node definitions             [x y z] (3 * 3)
%   dNodedx    Node definitions derivatives         (SIZE(Node) * nVar)
%   t          Transformation matrix                (3 * 3)
%   dtdx       Transformation matrix derivatives    (3 * 3 * nVar)
%
%   See also KE_BEAM, TRANS_TRUSS.

% David Dooms, Wouter Dillen
% March 2008, April 2017

if nargin<2, dNodedx = []; end
nVar = 0;
if nargout>1 && ~isempty(dNodedx)
    nVar = size(dNodedx,3);
end
if nVar==0 || isempty(dNodedx), dNodedx = zeros([size(Node),nVar]); end


Nx=Node(2,:)-Node(1,:);        % vector along the local x-axis
Nx=Nx/norm(Nx);                % normalized vector along the local x-axis

V13=Node(3,:)-Node(1,:);       % vector between nodes 1 and 3

Nz=cross(Nx,V13);              % vector along the local z-axis
if Nz==0
    error('The three nodes of the beam element are collinear')
end
Nz=Nz/norm(Nz);                % normalized vector along the local z-axis
Ny=cross(Nz,Nx);               % normalized vector along the local y-axis

t=[Nx; Ny; Nz];


dtdx=zeros(3,3,nVar);
for n=1:nVar
    if any(reshape(dNodedx(:,:,n),[],1)~=0)
        Nx=Node(2,:)-Node(1,:);
        dNxdx=dNodedx(2,:,n)-dNodedx(1,:,n);
        dNxdx=normder(Nx,dNxdx);
        Nx=Nx/norm(Nx);
        
        V13=Node(3,:)-Node(1,:);
        dV13dx=dNodedx(3,:,n)-dNodedx(1,:,n);
        
        Nz=cross(Nx,V13);        
        dNzdx=cross(dNxdx,V13)+cross(Nx,dV13dx);
        dNzdx=normder(Nz,dNzdx);       
        Nz=Nz/norm(Nz);
        dNydx=cross(dNzdx,Nx)+cross(Nz,dNxdx);
        
        dtdx(:,:,n)=[dNxdx; dNydx; dNzdx];
    end
end

end


function der=normder(v,dv)
    % returns the derivative of (row) vector normalization d(v/norm(v))/dx
    der=(dv*(v*v.')-v*(dv*v.'))/(norm(v)^3);
end
    
