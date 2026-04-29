function [Ke,Me,dKedx] = ke_truss(Node,Section,Material,Options,dNodedx,dSectiondx)

%KE_TRUSS   Truss element stiffness and mass matrix in global coordinate system.
%
%   [Ke,Me] = KE_TRUSS(Node,Section,Material,Options) returns the element
%   stiffness and mass matrix in the global coordinate system 
%   for a two node truss element (isotropic material).
%
%   [Ke,~,dKedx] = KE_TRUSS(Node,Section,Material,Options,dNodedx,dSectiondx) 
%   returns the element stiffness matrix in the global coordinate system
%   for a two node truss element (isotropic material), and additionally
%   computes the derivatives of the stiffness matrix with respect to the
%   design variables x. The derivatives of the mass matrix have not yet
%   been implemented.
%
%   Node       Node definitions           [x1 y1 z1; x2 y2 z2] (2 * 3)
%   Section    Section definition         [A]
%   Material   Material definition        [E nu rho]
%   Options    Struct containing element options. Fields: 
%      .lumped Construct lumped mass matrix:
%              {true | false (default)}
%   dNodedx    Node definitions derivatives     (SIZE(Node) * nVar)
%   dSectiondx Section definitions derivatives  (SIZE(Section) * nVar)
%   Ke         Element stiffness matrix         (6 * 6)
%   Me         Element mass matrix              (6 * 6)
%   dKedx      Element stiffness matrix         (CELL(nVar,1))
%
%   See also KELCS_TRUSS, TRANS_TRUSS, ASMKM, KE_BEAM.

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
if ~ all(isfinite(Node(1:2,1:3)))
    error('Not all the nodes exist.')
end

% Element length
L=norm(Node(2,:)-Node(1,:)); 
dLdx=(permute(dNodedx(2,:,:)-dNodedx(1,:,:),[3 2 1])*(Node(2,:)-Node(1,:)).')/L;

% Material
E=Material(1);

% Section
A=Section(1);
dAdx=dSectiondx(:,1,:);

if nargout==2           % stiffness and mass
    
    t=trans_truss(Node);
    T=blkdiag(t,t);
    %
    rho=Material(3);
    [KeLCS,MeLCS]=kelcs_truss(L,A,E,rho,Options);
    Me=full(MeLCS);
    Ke=T.'*KeLCS*T;
    
else                    % only stiffness
    
    [t,dtdx]=trans_truss(Node,dNodedx);
    T=blkdiag(t,t);
    dTdx = [dtdx  zeros(size(dtdx));
            zeros(size(dtdx)) dtdx];
    %
    [KeLCS,~,dKeLCSdx]=kelcs_truss(L,A,E,[],{},dLdx,dAdx);
    Ke=T.'*KeLCS*T;
    for n=1:nVar
        if any(reshape(dtdx(:,:,n),[],1)~=0)
            dKedx{n} = dTdx(:,:,n).'*(KeLCS*T) + T.'*(dKeLCSdx{n}*T + KeLCS*dTdx(:,:,n));
        else
            dKedx{n} = T.'*dKeLCSdx{n}*T;
        end
    end 
    
end
