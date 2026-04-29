function [u,t]=cdiff(M,C,K,dt,p,u0,u1)

%CDIFF   Direct time integration for dynamic systems - central diff. method.
%   [u,t] = CDIFF(M,C,K,dt,p,u0,u1) applies the central difference method for
%   the calculation of the nodal displacements u of the dynamic system with
%   the system matrices M, C and K due to the excitation p.
%
%   M    Mass matrix (nDof * nDof)
%   C    Damping matrix (nDof * nDof)
%   K    Stiffness matrix (nDof * nDof)
%   dt   Time step of the integration scheme (1 * 1).  Should be small enough
%        to ensure the stability and the precision of the integration scheme.
%   p    Excitation (nDof * N).  p(:,k) corresponds to time point t(k).
%   u0   Displacements at time point t(1)-dt (nDof * 1).  Defaults to zero.
%   u1   Displacements at time point t(1) (nDof * 1).  Defaults to zero.
%   u    Displacements (nDof * N).  u(:,k) corresponds to time point t(k).
%   t    Time axis (1 * N), defined as t = [0:N-1] * dt.

% Mattias Schevenels
% March 2004

% INITIALISATION
nDof=size(M,1);
N=size(p,2);
t=[0:N-1]*dt;

% DEFAULT INPUT ARGUMENTS
if nargin<6, u0=zeros(nDof,1); end;
if nargin<7, u1=zeros(nDof,1); end;

% CHECK INPUT ARGUMENTS
error(argdimchk( ...
  M, {'nDof','nDof'}, ...
  C, {'nDof','nDof'}, ...
  K, {'nDof','nDof'}, ...
  dt,{ 1    , 1    }, ...
  p, {'nDof','N'   }, ...
  u0,{'nDof', 1    }, ...
  u1,{'nDof', 1    }));

% DIRECT INTEGRATION: CENTRAL DIFFERENCE METHOD
% (Course text, p. 7.19-7.20)
L=inv(M/dt^2+C/2/dt);
R1=2*M/dt^2-K;
R2=-M/dt^2+C/2/dt;
u=zeros(nDof,N);
u(:,1)=u1;

for k=2:N
  if k==2
    R=p(:,k-1)+R1*u(:,k-1)+R2*u0;
    u(:,k)=L*R;
  else
    R=p(:,k-1)+R1*u(:,k-1)+R2*u(:,k-2);
    u(:,k)=L*R;
  end
end
