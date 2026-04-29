function [C,C_lambda,C_mu]=cmat_isotropic(problem,Section,Material)

%CMAT_ISOTROPIC  Constitutive matrix for isotropic materials.
%
%   [C,C_lambda,C_mu]==cmat_isotropic(problem,Section,Material)
%   computes the constitutive matrix for isotropic materials.
%
%	problem		
%	Section		Section definition
%	Material	Material definition
%	C			Constitutive matrix (nStress * nStrain)
%   C_lambda    Contribution of lambda to C (nStress * nStrain)
%   C_mu        Contribution of mu to C (nStress * nStrain)

% Jef Wambacq
% 2017

E=Material(1);
nu=Material(2);

lambda=nu*E/((1+nu)*(1-2*nu));
mu=E/(2*(1+nu));

switch lower(problem)
    case '2dstress'
        if isempty(Section)
            h=1;
        else
            h=Section(1);
        end
        
        % Cannot be split in a lambda and mu contribution
        C=E*h/(1-nu^2)*[1  nu 0;
                        nu 1  0;
                        0  0  (1-nu)/2];
    case '2dstrain'
        C_lambda=lambda*[1 1 0;
                         1 1 0;
                         0 0 0];
        C_mu=mu*[2 0 0;
                 0 2 0;
                 0 0 1];
        
        C=C_lambda+C_mu;
	case 'axisym'
        C_lambda=lambda*[1 1 1 0;
                         1 1 1 0;
                         1 1 1 0;
                         0 0 0 0];
        C_mu=mu*[2 0 0 0;
                 0 2 0 0;
                 0 0 2 0;
                 0 0 0 1];
                     
        C=C_lambda+C_mu;
    otherwise
        C_lambda=lambda*[1 1 1 0 0 0;
                         1 1 1 0 0 0;
                         1 1 1 0 0 0;
                         0 0 0 0 0 0;
                         0 0 0 0 0 0;
                         0 0 0 0 0 0];
        C_mu=mu*[2 0 0 0 0 0;
                 0 2 0 0 0 0;
                 0 0 2 0 0 0;
                 0 0 0 1 0 0;
                 0 0 0 0 1 0;
                 0 0 0 0 0 1];
                     
        C=C_lambda+C_mu;
end
end
