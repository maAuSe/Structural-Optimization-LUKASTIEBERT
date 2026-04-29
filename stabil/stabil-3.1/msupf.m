function [X,H]=msupf(omega,xi,Omega,Pm);

%MSUPF   Modal superposition in the frequency domain.
%
%   [X,H] = MSUPF(omega,xi,Omega,Pm) calculates the modal transfer functions H
%   and the modal displacements x(t) = X * EXP(i * Omega * t) of a dynamic
%   system with the eigenfrequencies omega and the modal damping ratios xi due
%   to the modal excitation pm(t) = Pm * EXP(i * Omega * t).
%
%   omega   Eigenfrequencies [rad/s] (nMode * 1)
%   xi      Modal damping ratios [ - ] (nMode * 1) or (1 * 1), constant modal
%           damping is assumed in the (1 * 1) case.
%   Omega   Excitation frequencies [rad/s] (1 * N).
%   Pm      Complex amplitude of the modal excitation (nMode * N).
%   X       Complex amplitude of the modal displacements (nMode * N).
%   H       Modal transfer functions (nMode * N).

% Mattias Schevenels
% March 2004

% INITIALISATION
nMode=length(omega); % Number of modes to take into account
N=size(Pm,2);        % Number of excitation frequencies

% DEFAULT INPUT ARGUMENTS
if length(xi)==1, xi=repmat(xi,nMode,1); end; % Constant modal damping

% CHECK INPUT ARGUMENTS
error(argdimchk( ...
  omega,{'nMode', 1     }, ...
  xi,   {'nMode', 1     }, ...
  Omega,{ 1     ,'N'    }, ...
  Pm,   {'nMode','N'    }));

% TRANSFER FUNCTION (nMode * N)
H=1./( repmat(-Omega.^2,nMode,1) + 2*i*xi.*omega*Omega + repmat(omega.^2,1,N) );

% MODAL DISPLACEMENTS (nMode * N)
X=H.*Pm;

