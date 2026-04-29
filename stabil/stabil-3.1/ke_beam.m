function [Ke,Me,dKedx] = ke_beam(Node,Section,Material,Options,dNodedx,dSectiondx)

%KE_BEAM   Beam element stiffness and mass matrix in global coordinate system.
%
%   [Ke,Me] = KE_BEAM(Node,Section,Material,Options) returns the element
%   stiffness and mass matrix in the global coordinate system 
%   for a two node beam element (isotropic material).
%
%   [Ke,~,dKedx] = KE_BEAM(Node,Section,Material,Options,dNodedx,dSectiondx) 
%   returns the element stiffness matrix in the global coordinate system
%   for a two node beam element (isotropic material), and additionally
%   computes the derivatives of the stiffness matrix with respect to the
%   design variables x. The derivatives of the mass matrix have not yet
%   been implemented.
%
%   Node       Node definitions           [x y z] (3 * 3)
%   Section    Section definition         [A ky kz Ixx Iyy Izz]
%   Material   Material definition        [E nu rho]
%   Options    Struct containing element options. Fields: 
%      .lumped            Construct lumped mass matrix:
%                         {true | false (default)}
%      .rotationalInertia Include rotational inertia:
%                         {true (default)| false}
%   dNodedx    Node definitions derivatives         (SIZE(Node) * nVar)
%   dSectiondx Section definitions derivatives      (SIZE(Section) * nVar)
%   Ke         Element stiffness matrix             (12 * 12)
%   Me         Element mass matrix                  (12 * 12)
%   dKedx      Element stiffness matrix derivatives (CELL(nVar,1))
%
%   See also KELCS_BEAM, TRANS_BEAM, ASMKM, KE_TRUSS.

% David Dooms, Wouter Dillen
% March 2008, April 2017

if nargin<4, Options = {}; end
if nargin<5, dNodedx = []; end
if nargin<6, dSectiondx = []; end

nVar = 0;
if nargout>2 && (~isempty(dNodedx) || ~isempty(dSectiondx))
    nVar = max(size(dNodedx,3),size(dSectiondx,3));
end

if nVar==0 || isempty(dNodedx), dNodedx = zeros([size(Node),nVar]); end
if nVar==0 || isempty(dSectiondx), dSectiondx = zeros([size(Section),nVar]); end

if nargout>1, Me = []; end
if nargout>2, dKedx = cell(nVar,1); end


% Check nodes
if ~ all(isfinite(Node(1:3,1:3)))
    error('Not all the nodes exist.')
end

% Element length
L=norm(Node(2,:)-Node(1,:)); 
dLdx=(permute(dNodedx(2,:,:)-dNodedx(1,:,:),[3 2 1])*(Node(2,:)-Node(1,:)).')/L;

% Material
E=Material(1);
nu=Material(2);

% Section
A=Section(1);
ky=Section(2);
kz=Section(3);
Ixx=Section(4);
Iyy=Section(5);
Izz=Section(6);
dAdx=dSectiondx(:,1,:);
dkydx=dSectiondx(:,2,:);
dkzdx=dSectiondx(:,3,:);
dIxxdx=dSectiondx(:,4,:);
dIyydx=dSectiondx(:,5,:);
dIzzdx=dSectiondx(:,6,:);

if nargout==2           % stiffness and mass
    
    t=trans_beam(Node); 
    T=blkdiag(t,t,t,t);
    %
    rho=Material(3);
    [KeLCS,MeLCS]=kelcs_beam(L,A,ky,kz,Ixx,Iyy,Izz,E,nu,rho,Options);
    Me=T.'*MeLCS*T;
    Ke=T.'*KeLCS*T;
    
else                    % only stiffness
    
    [t,dtdx]=trans_beam(Node,dNodedx); 
    T=blkdiag(t,t,t,t);
    dTdx = [dtdx                zeros(size(dtdx)) 	zeros(size(dtdx))   zeros(size(dtdx));
            zeros(size(dtdx))   dtdx                zeros(size(dtdx))   zeros(size(dtdx));
            zeros(size(dtdx))   zeros(size(dtdx))   dtdx                zeros(size(dtdx));
            zeros(size(dtdx))   zeros(size(dtdx))   zeros(size(dtdx))   dtdx];
    %
    [KeLCS,~,dKeLCSdx]=kelcs_beam(L,A,ky,kz,Ixx,Iyy,Izz,E,nu,[],{},dLdx,dAdx,dkydx,dkzdx,dIxxdx,dIyydx,dIzzdx);
    Ke=T.'*KeLCS*T;
    for n=1:nVar
        if any(reshape(dtdx(:,:,n),[],1)~=0)
            dKedx{n} = dTdx(:,:,n).'*(KeLCS*T) + T.'*(dKeLCSdx{n}*T + KeLCS*dTdx(:,:,n));
        else
            dKedx{n} = T.'*dKeLCSdx{n}*T;
        end
    end
    
end
