function [KeLCS,MeLCS,dKeLCSdx] = kelcs_truss(L,A,E,rho,Options,dLdx,dAdx)

%KELCS_TRUSS   Truss element stiffness and mass matrix in local coordinate system.
%
%   [KeLCS,MeLCS] = KELCS_TRUSS(L,A,E,rho,Options) 
%   [KeLCS,MeLCS] = KELCS_TRUSS(L,A,E,rho) 
%    KeLCS        = KELCS_TRUSS(L,A,E)
%   returns the element stiffness and mass matrix in the local coordinate system  
%   for a two node truss element (isotropic material).
%
%   [KeLCS,~,dKeLCSdx] = KELCS_TRUSS(L,A,E,rho,Options,dLdx,dAdx) 
%   returns the element stiffness matrix in the local coordinate system
%   for a two node truss element (isotropic material), and additionally
%   computes the derivatives of the stiffness matrix with respect to the
%   design variables x. The derivatives of the mass matrix have not yet
%   been implemented.
%
%   L          Truss length 
%   A          Truss cross section area
%   E          Young's modulus
%   rho        Mass density
%   Options    Options for the mass matrix: {'lumped'}
%   dLdx       Truss length derivatives                 (1 * nVar)
%   dAdx       Truss cross section area derivatives     (1 * nVar)
%   KeLCS      Element stiffness matrix                 (6 * 6)
%   MeLCS      Element mass matrix                      (6 * 6)
%   dKeLCSdx   Element stiffness matrix derivatives     (CELL(nVar,1))
%
%   See also KE_TRUSS, KELCS_BEAM.

% David Dooms, Wouter Dillen
% March 2008, April 2017

if nargin<5, Options = {}; end
if nargin<6, dLdx = []; end
if nargin<7, dAdx = []; end

nVar = 0;
if nargout>2 && (~isempty(dLdx) || ~isempty(dAdx))
    nVar = max(length(dLdx),length(dAdx));
end

if nVar==0 || isempty(dLdx), dLdx = zeros(1,1,nVar); end
if nVar==0 || isempty(dAdx), dAdx = zeros(1,1,nVar); end

KeLCS=sparse([],[],[],6,6,4);
if nargout>1, MeLCS=sparse([]); end
if nargout>2
    dKeLCSdx=cell(nVar,1);
    dKeLCSdx_temp1=sparse([],[],[],6,6,4);   % section properties derivatives
    dKeLCSdx_temp2=sparse([],[],[],6,6,4);   % geometry derivatives
end


%% STIFFNESS MATRIX %%

KeLCS(1,1) = E*A/L;
KeLCS(4,4)=KeLCS(1,1);
KeLCS(4,1)=-KeLCS(1,1);
% symmetry
KeLCS(1,4)=KeLCS(4,1);

if nVar>0
    
    % TEMP1 - section properties derivatives
    dKeLCSdx_temp1(1,1) = E/L;
    dKeLCSdx_temp1(4,4)=dKeLCSdx_temp1(1,1);
    dKeLCSdx_temp1(4,1)=-dKeLCSdx_temp1(1,1);
    % symmetry
    dKeLCSdx_temp1(1,4)=dKeLCSdx_temp1(4,1);
    
    % TEMP2 - geometry derivatives
    dKeLCSdx_temp2(1,1) = -E*A/L^2;
    dKeLCSdx_temp2(4,4)=dKeLCSdx_temp2(1,1);
    dKeLCSdx_temp2(4,1)=-dKeLCSdx_temp2(1,1);
    % symmetry
    dKeLCSdx_temp2(1,4)=dKeLCSdx_temp2(4,1);
    
end

for n=1:nVar
    dKeLCSdx{n} = dKeLCSdx_temp1*dAdx(n) + dKeLCSdx_temp2*dLdx(n);
end

%% MASS MATRIX %%

if nargout==2
    
    if ~isfield(Options,'lumped')
        Options.lumped=false;
    end

    if ~Options.lumped      % consistent mass matrix

        MeLCS=sparse([],[],[],6,6,12);

        MeLCS(1,1)=rho*A*L/3;
        MeLCS(2,2)=MeLCS(1,1);
        MeLCS(3,3)=MeLCS(1,1);
        MeLCS(4,4)=MeLCS(1,1);
        MeLCS(5,5)=MeLCS(1,1);
        MeLCS(6,6)=MeLCS(1,1);

        MeLCS(4,1)=rho*A*L/6;
        MeLCS(5,2)=MeLCS(4,1);
        MeLCS(6,3)=MeLCS(4,1);

        % symmetry
        MeLCS(1,4)=MeLCS(4,1);
        MeLCS(2,5)=MeLCS(5,2);
        MeLCS(3,6)=MeLCS(6,3);

    else                                    % lumped mass matrix
        
        MeLCS=sparse([],[],[],6,6,6);

        MeLCS(1,1)=rho*A*L/2;
        MeLCS(2,2)=MeLCS(1,1);
        MeLCS(3,3)=MeLCS(1,1);
        MeLCS(4,4)=MeLCS(1,1);
        MeLCS(5,5)=MeLCS(1,1);
        MeLCS(6,6)=MeLCS(1,1);
        
    end

end
