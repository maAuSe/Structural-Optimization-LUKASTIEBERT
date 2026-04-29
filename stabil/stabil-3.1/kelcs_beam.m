function [KeLCS,MeLCS,dKeLCSdx] = kelcs_beam(L,A,ky,kz,Ixx,Iyy,Izz,E,nu,rho,Options,dLdx,dAdx,dkydx,dkzdx,dIxxdx,dIyydx,dIzzdx)

%KELCS_BEAM   Beam element stiffness and mass matrix in local coordinate system.
%
%   [KeLCS,MeLCS] = KELCS_BEAM(L,A,ky,kz,Ixx,Iyy,Izz,E,nu,rho,Options) 
%   [KeLCS,MeLCS] = KELCS_BEAM(L,A,ky,kz,Ixx,Iyy,Izz,E,nu,rho) 
%    KeLCS        = KELCS_BEAM(L,A,ky,kz,Ixx,Iyy,Izz,E,nu)
%   returns the element stiffness and mass matrix in the local coordinate system 
%   for a two node beam element (isotropic material).
%
%   [KeLCS,~,dKeLCSdx] = KELCS_BEAM(L,A,ky,kz,Ixx,Iyy,Izz,E,nu,rho,Options,dLdx,dAdx,...
%                                                            dkydx,dkzdx,dIxxdx,dIyydx,dIzzdx) 
%   returns the element stiffness matrix in the local coordinate system
%   for a two node beam element (isotropic material), and additionally
%   computes the derivatives of the stiffness matrix with respect to the
%   design variables x. The derivatives of the mass matrix have not yet
%   been implemented.
%
%   L          Beam length 
%   A          Beam cross section area
%   ky         Shear deflection factor A_sy = ky * A
%   kz         Shear deflection factor A_sz = kz * A
%   Ixx        Moment of inertia
%   Iyy        Moment of inertia
%   Izz        Moment of inertia
%   E          Young's modulus
%   nu         Poisson coefficient
%   rho        Mass density
%   Options    Options for the mass matrix: {'lumped'}, {'norotaroryinertia'}
%   dLdx       Beam length derivatives              (1 * nVar)
%   dAdx       Beam cross section area derivatives  (1 * nVar)
%   dkydx      Shear deflection factor derivatives  (1 * nVar)
%   dkzdx      Shear deflection factor derivatives  (1 * nVar)
%   dIxxdx     Moment of inertia derivatives        (1 * nVar)
%   dIyydx     Moment of inertia derivatives        (1 * nVar)
%   dIzzdx     Moment of inertia derivatives        (1 * nVar)
%   KeLCS      Element stiffness matrix             (12 * 12)
%   MeLCS      Element mass matrix                  (12 * 12)
%   dKeLCSdx   Element stiffness matrix derivatives (CELL(nVar,1))
%
%   See also KE_BEAM, KELCS_TRUSS.

% David Dooms, Wouter Dillen
% March 2008, April 2017

%%% CONTROLES %%% alles groter dan nul, nu kleiner dan 0.5
% bij ky en kz wordt heel kleine waarde opgeteld om problemen te vermijden als ze gelijk zijn aan nul
% default for Ixx?
% added mass?

if nargin<11, Options = {}; end
if nargin<12, dLdx = []; end
if nargin<13, dAdx = []; end
if nargin<14, dkydx = []; end
if nargin<15, dkzdx = []; end
if nargin<16, dIxxdx = []; end
if nargin<17, dIyydx = []; end
if nargin<18, dIzzdx = []; end

nVar = 0;
if nargout>2 && (~isempty(dLdx) || ~isempty(dAdx) || ~isempty(dkydx) || ~isempty(dkzdx) || ~isempty(dIxxdx) || ~isempty(dIyydx) || ~isempty(dIzzdx))
    nVar = max([length(dLdx),length(dAdx),length(dkydx),length(dkzdx),length(dIxxdx),length(dIyydx),length(dIzzdx)]);
end

if nVar==0 || isempty(dLdx), dLdx = zeros(1,1,nVar); else, dLdx = permute(dLdx(:),[2 3 1]); end
if nVar==0 || isempty(dAdx), dAdx = zeros(1,1,nVar); else, dAdx = permute(dAdx(:),[2 3 1]); end
if ~isempty(dkydx) && any(dkydx(:))~=0, error('Sensitivities for the shear deflection factor (ky) have not been implemented yet.'); end
if ~isempty(dkzdx) && any(dkzdx(:))~=0, error('Sensitivities for the shear deflection factor (kz) have not been implemented yet.'); end
if nVar==0 || isempty(dIxxdx), dIxxdx = zeros(1,1,nVar); else, dIxxdx = permute(dIxxdx(:),[2 3 1]); end
if nVar==0 || isempty(dIyydx), dIyydx = zeros(1,1,nVar); else, dIyydx = permute(dIyydx(:),[2 3 1]); end
if nVar==0 || isempty(dIzzdx), dIzzdx = zeros(1,1,nVar); else, dIzzdx = permute(dIzzdx(:),[2 3 1]); end

KeLCS=sparse([],[],[],12,12,40);
if nargout>1, MeLCS=sparse([]); end
if nargout>2
    dKeLCSdx=cell(nVar,1);
    dKeLCSdx_temp1=zeros(12,12,nVar);   % section property derivatives
    dKeLCSdx_temp2=zeros(12,12);        % geometry derivatives
end


%% STIFFNESS MATRIX %%
% Przemieniecki, J. S., Theory of Matrix Structural Analysis, McGraw-Hill,
% New York (1968).

G=E/(2*(1+nu));

phi_y=(12*E*Izz)/(G*A*(ky+1e-250)*L^2);
phi_z=(12*E*Iyy)/(G*A*(kz+1e-250)*L^2);

% terms with A
KeLCS(1,1)=A*E/L;
KeLCS(7,7)=KeLCS(1,1);
KeLCS(7,1)=-KeLCS(1,1);

% terms with Ixx
KeLCS(4,4)=Ixx*G/L;
KeLCS(10,10)=KeLCS(4,4);
KeLCS(10,4)=-KeLCS(4,4);
                     
% terms with Izz
KeLCS(2,2)=(12*E*Izz)/(L^3*(1+phi_y));
KeLCS(8,8)=KeLCS(2,2);
KeLCS(8,2)=-KeLCS(2,2);

KeLCS(6,2)=(6*E*Izz)/(L^2*(1+phi_y));
KeLCS(12,2)=KeLCS(6,2);
KeLCS(8,6)=-KeLCS(6,2);
KeLCS(12,8)=-KeLCS(6,2);

KeLCS(6,6)=E*Izz/L*(4+phi_y)/(1+phi_y);
KeLCS(12,12)=KeLCS(6,6);

KeLCS(12,6)=E*Izz/L*(2-phi_y)/(1+phi_y);

% terms with Iyy
KeLCS(3,3)=(12*E*Iyy)/(L^3*(1+phi_z));
KeLCS(9,9)=KeLCS(3,3);
KeLCS(9,3)=-KeLCS(3,3);

KeLCS(9,5)=(6*E*Iyy)/(L^2*(1+phi_z));
KeLCS(11,9)=KeLCS(9,5);
KeLCS(5,3)=-KeLCS(9,5);
KeLCS(11,3)=-KeLCS(9,5);

KeLCS(5,5)=E*Iyy/L*(4+phi_z)/(1+phi_z);
KeLCS(11,11)=KeLCS(5,5);

KeLCS(11,5)=E*Iyy/L*(2-phi_z)/(1+phi_z);

% symmetry
KeLCS(1,7)=KeLCS(7,1);
KeLCS(4,10)=KeLCS(10,4);
KeLCS(2,8)=KeLCS(8,2);
KeLCS(2,6)=KeLCS(6,2);
KeLCS(2,12)=KeLCS(12,2);
KeLCS(6,8)=KeLCS(8,6);
KeLCS(8,12)=KeLCS(12,8);
KeLCS(6,12)=KeLCS(12,6);
KeLCS(3,9)=KeLCS(9,3);
KeLCS(5,9)=KeLCS(9,5);
KeLCS(9,11)=KeLCS(11,9);
KeLCS(3,5)=KeLCS(5,3);
KeLCS(3,11)=KeLCS(11,3);
KeLCS(5,11)=KeLCS(11,5);


if nVar>0
    
    if isfinite(ky)  &&  any([dAdx(:); dIzzdx(:); dLdx(:)]) > 0
        % dphi_ydx is not equal to zero
        error('Sensitivities for shear deformation in beam elements have not been implemented yet.')
    end
    if isfinite(kz)  &&  any([dAdx(:); dIyydx(:); dLdx(:)]) > 0
        % dphi_zdx is not equal to zero
        error('Sensitivities for shear deformation in beam elements have not been implemented yet.')
    end
        
    % TEMP1 - section properties derivatives
    if any([dAdx(:); dIxxdx(:); dIyydx(:); dIzzdx(:)]) > 0
                        
    % terms with A
    dKeLCSdx_temp1(1,1,:)=dAdx*E/L;
    dKeLCSdx_temp1(7,7,:)=dKeLCSdx_temp1(1,1,:);
    dKeLCSdx_temp1(7,1,:)=-dKeLCSdx_temp1(1,1,:);

    % terms with Ixx
    dKeLCSdx_temp1(4,4,:)=dIxxdx*G/L;
    dKeLCSdx_temp1(10,10,:)=dKeLCSdx_temp1(4,4,:);
    dKeLCSdx_temp1(10,4,:)=-dKeLCSdx_temp1(4,4,:);
                     
    % terms with Izz
    dKeLCSdx_temp1(2,2,:)=(12*E*dIzzdx)/(L^3*(1+phi_y));
    dKeLCSdx_temp1(8,8,:)=dKeLCSdx_temp1(2,2,:);
    dKeLCSdx_temp1(8,2,:)=-dKeLCSdx_temp1(2,2,:);

    dKeLCSdx_temp1(6,2,:)=(6*E*dIzzdx)/(L^2*(1+phi_y));
    dKeLCSdx_temp1(12,2,:)=dKeLCSdx_temp1(6,2,:);
    dKeLCSdx_temp1(8,6,:)=-dKeLCSdx_temp1(6,2,:);
    dKeLCSdx_temp1(12,8,:)=-dKeLCSdx_temp1(6,2,:);

    dKeLCSdx_temp1(6,6,:)=E*dIzzdx/L*(4+phi_y)/(1+phi_y);
    dKeLCSdx_temp1(12,12,:)=dKeLCSdx_temp1(6,6,:);

    dKeLCSdx_temp1(12,6,:)=E*dIzzdx/L*(2-phi_y)/(1+phi_y);

    % terms with Iyy
    dKeLCSdx_temp1(3,3,:)=(12*E*dIyydx)/(L^3*(1+phi_z));
    dKeLCSdx_temp1(9,9,:)=dKeLCSdx_temp1(3,3,:);
    dKeLCSdx_temp1(9,3,:)=-dKeLCSdx_temp1(3,3,:);

    dKeLCSdx_temp1(9,5,:)=(6*E*dIyydx)/(L^2*(1+phi_z));
    dKeLCSdx_temp1(11,9,:)=dKeLCSdx_temp1(9,5,:);
    dKeLCSdx_temp1(5,3,:)=-dKeLCSdx_temp1(9,5,:);
    dKeLCSdx_temp1(11,3,:)=-dKeLCSdx_temp1(9,5,:);

    dKeLCSdx_temp1(5,5,:)=E*dIyydx/L*(4+phi_z)/(1+phi_z);
    dKeLCSdx_temp1(11,11,:)=dKeLCSdx_temp1(5,5,:);

    dKeLCSdx_temp1(11,5,:)=E*dIyydx/L*(2-phi_z)/(1+phi_z);

    % symmetry
    dKeLCSdx_temp1(1,7,:)=dKeLCSdx_temp1(7,1,:);
    dKeLCSdx_temp1(4,10,:)=dKeLCSdx_temp1(10,4,:);
    dKeLCSdx_temp1(2,8,:)=dKeLCSdx_temp1(8,2,:);
    dKeLCSdx_temp1(2,6,:)=dKeLCSdx_temp1(6,2,:);
    dKeLCSdx_temp1(2,12,:)=dKeLCSdx_temp1(12,2,:);
    dKeLCSdx_temp1(6,8,:)=dKeLCSdx_temp1(8,6,:);
    dKeLCSdx_temp1(8,12,:)=dKeLCSdx_temp1(12,8,:);
    dKeLCSdx_temp1(6,12,:)=dKeLCSdx_temp1(12,6,:);
    dKeLCSdx_temp1(3,9,:)=dKeLCSdx_temp1(9,3,:);
    dKeLCSdx_temp1(5,9,:)=dKeLCSdx_temp1(9,5,:);
    dKeLCSdx_temp1(9,11,:)=dKeLCSdx_temp1(11,9,:);
    dKeLCSdx_temp1(3,5,:)=dKeLCSdx_temp1(5,3,:);
    dKeLCSdx_temp1(3,11,:)=dKeLCSdx_temp1(11,3,:);
    dKeLCSdx_temp1(5,11,:)=dKeLCSdx_temp1(11,5,:);
    end
    
    
    % TEMP2 - geometry derivatives
    % terms with A
    dKeLCSdx_temp2(1,1)=A*E/L^2;
    dKeLCSdx_temp2(7,7)=dKeLCSdx_temp2(1,1);
    dKeLCSdx_temp2(7,1)=-dKeLCSdx_temp2(1,1);

    % terms with Ixx
    dKeLCSdx_temp2(4,4)=Ixx*G/L^2;
    dKeLCSdx_temp2(10,10)=dKeLCSdx_temp2(4,4);
    dKeLCSdx_temp2(10,4)=-dKeLCSdx_temp2(4,4);
                     
    % terms with Izz
    dKeLCSdx_temp2(2,2)=12*E*Izz/L^4*(3/(1+phi_y)-2*phi_y/(1+phi_y)^2);
    dKeLCSdx_temp2(8,8)=dKeLCSdx_temp2(2,2);
    dKeLCSdx_temp2(8,2)=-dKeLCSdx_temp2(2,2);

    dKeLCSdx_temp2(6,2)=6*E*Izz/L^3*(2/(1+phi_y)-2*phi_y/(1+phi_y)^2);
    dKeLCSdx_temp2(12,2)=dKeLCSdx_temp2(6,2);
    dKeLCSdx_temp2(8,6)=-dKeLCSdx_temp2(6,2);
    dKeLCSdx_temp2(12,8)=-dKeLCSdx_temp2(6,2);

    dKeLCSdx_temp2(6,6)=E*Izz/L^2*((4+phi_y)/(1+phi_y)-6*phi_y/(1+phi_y)^2);
    dKeLCSdx_temp2(12,12)=dKeLCSdx_temp2(6,6);

    dKeLCSdx_temp2(12,6)=E*Izz/L^2*((2-phi_y)/(1+phi_y)-6*phi_y/(1+phi_y)^2);

    % terms with Iyy
    dKeLCSdx_temp2(3,3)=12*E*Iyy/L^4*(3/(1+phi_z)-2*phi_z/(1+phi_z)^2);
    dKeLCSdx_temp2(9,9)=dKeLCSdx_temp2(3,3);
    dKeLCSdx_temp2(9,3)=-dKeLCSdx_temp2(3,3);

    dKeLCSdx_temp2(9,5)=6*E*Iyy/L^3*(2/(1+phi_z)-2*phi_z/(1+phi_z)^2);
    dKeLCSdx_temp2(11,9)=dKeLCSdx_temp2(9,5);
    dKeLCSdx_temp2(5,3)=-dKeLCSdx_temp2(9,5);
    dKeLCSdx_temp2(11,3)=-dKeLCSdx_temp2(9,5);

    dKeLCSdx_temp2(5,5)=E*Iyy/L^2*((4+phi_z)/(1+phi_z)-6*phi_z/(1+phi_z)^2);
    dKeLCSdx_temp2(11,11)=dKeLCSdx_temp2(5,5);

    dKeLCSdx_temp2(11,5)=E*Iyy/L^2*((2-phi_z)/(1+phi_z)-6*phi_z/(1+phi_z)^2);

    % symmetry
    dKeLCSdx_temp2(1,7)=dKeLCSdx_temp2(7,1);
    dKeLCSdx_temp2(4,10)=dKeLCSdx_temp2(10,4);
    dKeLCSdx_temp2(2,8)=dKeLCSdx_temp2(8,2);
    dKeLCSdx_temp2(2,6)=dKeLCSdx_temp2(6,2);
    dKeLCSdx_temp2(2,12)=dKeLCSdx_temp2(12,2);
    dKeLCSdx_temp2(6,8)=dKeLCSdx_temp2(8,6);
    dKeLCSdx_temp2(8,12)=dKeLCSdx_temp2(12,8);
    dKeLCSdx_temp2(6,12)=dKeLCSdx_temp2(12,6);
    dKeLCSdx_temp2(3,9)=dKeLCSdx_temp2(9,3);
    dKeLCSdx_temp2(5,9)=dKeLCSdx_temp2(9,5);
    dKeLCSdx_temp2(9,11)=dKeLCSdx_temp2(11,9);
    dKeLCSdx_temp2(3,5)=dKeLCSdx_temp2(5,3);
    dKeLCSdx_temp2(3,11)=dKeLCSdx_temp2(11,3);
    dKeLCSdx_temp2(5,11)=dKeLCSdx_temp2(11,5);
        
end

for n=1:nVar
    dKeLCSdx{n} = sparse(dKeLCSdx_temp1(:,:,n)+dKeLCSdx_temp2*(-dLdx(n)));
end

%% MASS MATRIX %%

if nargout==2 
    
    if ~isfield(Options,'lumped')
        Options.lumped=false;
    end
    if ~isfield(Options,'rotationalInertia')
        Options.rotationalInertia=true;
    end

    if ~Options.lumped     % consistent mass matrix
    
        MeLCS=sparse([],[],[],12,12,40);

        % Yokoyama, T., "Vibrations of a Hanging Timoshenko Beam Under Gravity", Journal
        % of Sound and Vibration, Vol. 141, No. 2, pp. 245-258 (1990).

        if Options.rotationalInertia
            ry=sqrt(Iyy/A);
            rz=sqrt(Izz/A);
        else
            ry=0;
            rz=0;
        end

        % axial terms 
        MeLCS(1,1)=rho*A*L/3;
        MeLCS(7,7)=MeLCS(1,1);
        MeLCS(7,1)=rho*A*L/6;

        % torsional terms
        MeLCS(4,4)=rho*L*(Iyy+Izz)/3;
        MeLCS(10,10)=MeLCS(4,4);
        MeLCS(10,4)=rho*L*(Iyy+Izz)/6;

        % terms with Izz
        MeLCS(2,2)=rho*A*L/(1+phi_y)^2*(13/35+7/10*phi_y+1/3*phi_y^2+6/5*(rz/L)^2);
        MeLCS(8,8)=MeLCS(2,2);
        MeLCS(8,2)=rho*A*L/(1+phi_y)^2*(9/70+3/10*phi_y+1/6*phi_y^2-6/5*(rz/L)^2);

        MeLCS(6,2)=rho*A*L/(1+phi_y)^2*(11/210+11/120*phi_y+1/24*phi_y^2+(1/10-1/2*phi_y)*(rz/L)^2)*L;
        MeLCS(12,8)=-MeLCS(6,2); %controleren!!!!!
        MeLCS(8,6)=rho*A*L/(1+phi_y)^2*(13/420+3/40*phi_y+1/24*phi_y^2-(1/10-1/2*phi_y)*(rz/L)^2)*L;
        MeLCS(12,2)=-MeLCS(8,6);

        MeLCS(6,6)=rho*A*L/(1+phi_y)^2*(1/105+1/60*phi_y+1/120*phi_y^2+(2/15+1/6*phi_y+1/3*phi_y^2)*(rz/L)^2)*L^2;
        MeLCS(12,12)=MeLCS(6,6);

        MeLCS(12,6)=-rho*A*L/(1+phi_y)^2*(1/140+1/60*phi_y+1/120*phi_y^2+(1/30+1/6*phi_y-1/6*phi_y^2)*(rz/L)^2)*L^2;

        % terms with Iyy
        MeLCS(3,3)=rho*A*L/(1+phi_z)^2*(13/35+7/10*phi_z+1/3*phi_z^2+6/5*(ry/L)^2);
        MeLCS(9,9)=MeLCS(3,3);
        MeLCS(9,3)=rho*A*L/(1+phi_z)^2*(9/70+3/10*phi_z+1/6*phi_z^2-6/5*(ry/L)^2);

        MeLCS(5,3)=-rho*A*L/(1+phi_z)^2*(11/210+11/120*phi_z+1/24*phi_z^2+(1/10-1/2*phi_z)*(ry/L)^2)*L;
        MeLCS(11,9)=-MeLCS(5,3); %controleren!!!!!
        MeLCS(9,5)=-rho*A*L/(1+phi_z)^2*(13/420+3/40*phi_z+1/24*phi_z^2-(1/10-1/2*phi_z)*(ry/L)^2)*L;
        MeLCS(11,3)=-MeLCS(9,5);

        MeLCS(5,5)=rho*A*L/(1+phi_z)^2*(1/105+1/60*phi_z+1/120*phi_z^2+(2/15+1/6*phi_z+1/3*phi_z^2)*(ry/L)^2)*L^2;
        MeLCS(11,11)=MeLCS(5,5);

        MeLCS(11,5)=-rho*A*L/(1+phi_z)^2*(1/140+1/60*phi_z+1/120*phi_z^2+(1/30+1/6*phi_z-1/6*phi_z^2)*(ry/L)^2)*L^2;

        % symmetry
        MeLCS(1,7)=MeLCS(7,1);
        MeLCS(4,10)=MeLCS(10,4);
        MeLCS(2,8)=MeLCS(8,2);
        MeLCS(2,6)=MeLCS(6,2);
        MeLCS(2,12)=MeLCS(12,2);
        MeLCS(6,8)=MeLCS(8,6);
        MeLCS(8,12)=MeLCS(12,8);
        MeLCS(6,12)=MeLCS(12,6);
        MeLCS(3,9)=MeLCS(9,3);
        MeLCS(5,9)=MeLCS(9,5);
        MeLCS(9,11)=MeLCS(11,9);
        MeLCS(3,5)=MeLCS(5,3);
        MeLCS(3,11)=MeLCS(11,3);
        MeLCS(5,11)=MeLCS(11,5);

    else   % lumped 

        MeLCS=sparse([],[],[],12,12,12);

        MeLCS(1,1)=rho*A*L/2;
        MeLCS(2,2)=MeLCS(1,1);
        MeLCS(3,3)=MeLCS(1,1);
        MeLCS(7,7)=MeLCS(1,1);
        MeLCS(8,8)=MeLCS(1,1);
        MeLCS(9,9)=MeLCS(1,1);

        MeLCS(4,4)=rho*Ixx*L/2; % controleren!!!!!
        MeLCS(10,10)=MeLCS(4,4);

        MeLCS(5,5)=rho*A*L^3/24; % controleren!!!!! = rho*L/2*(Iyy/A+L^2/12); p594 Wunderlich
        MeLCS(6,6)=MeLCS(5,5);
        MeLCS(11,11)=MeLCS(5,5);
        MeLCS(12,12)=MeLCS(5,5);

    end

end

