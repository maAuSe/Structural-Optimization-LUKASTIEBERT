function C = cmat_shell8(MatType,Material,k)

%CMAT_SHELL8   Constitutive matrix for shell8 element
%
%   C = cmat_shell8(MatType,Material,k)
%   C = cmat_shell8(MatType,Material)
%   returns the constitutive matrix for shell8 elements.
%
%   MatType    Material type: 'isotropic' or 'orthotropic'
%   Material   Material definition
%              isotropic:   [E nu rho]
%              orthotropic: [Exx Eyy nuxy muxy muyz muzx theta rho]
%   k         Geometric coefficient for non-uniform shear stress (default = 1.2)
%   C         Constitutive matrix (5 * 5)
%
%   See also KE_SHELL8

% Miche Jansen
% 2012


if nargin<3
    k = 1.2;
end

switch lower(MatType)
    case 'orthotropic'
        % orthotropic: [Exx,Eyy,nuxy,muxy,muyz,muzx,theta,rho]
        Exx = Material(1,1);
        Eyy = Material(1,2);
        nuxy = Material(1,3);
        muxy = Material(1,4);
        muyz = Material(1,5);
        muzx = Material(1,6);
        theta = Material(1,7);
        
        B = Exx/(Exx-Eyy*nuxy^2);
        C = [Exx*B,Eyy*nuxy*B,0,0,0;
            B*Eyy*nuxy,Eyy*B,0,0,0;
            0,0,muxy,0,0;
            0,0,0,muyz/k,0;
            0,0,0,0,muzx/k];
        % transformatie matrix orthotroop assenkruis <-> lokaal
        if theta ~=0
            A = [cos(theta) sin(theta) 0;
                -sin(theta) cos(theta) 0;
                0 0 1];
            
            % transformatie matrix uit Cook p274 en p194 Zienkiewicz deel1
            T = vtrans_solid(A,'strain');
            T = T([1,2,4,5,6],[1,2,4,5,6]);
            % materiaalmatrix in globale assenkruis (p198 Zienkiewicz deel1)
            C = T.'*C*T;
        end
    case 'isotropic'
        % isotropic: material [E nu rho]
        E=Material(1,1);
        nu=Material(1,2);
        
        C = E/(1-nu^2)*[1 nu 0 0 0;         % is k = k uit cursus en niet kappa uit Zienkiewicz
            nu 1 0 0 0;
            0 0 (1-nu)/2 0 0;
            0 0 0 (1-nu)/(2*k) 0;
            0 0 0 0 (1-nu)/(2*k)];
end
end
