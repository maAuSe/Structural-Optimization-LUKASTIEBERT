function [u,v,a,t]=newmark(M,C,K,dt,p,varargin)

%NEWMARK   Direct time integration for dynamic systems - Newmark method
%   [u,v,a,t] = NEWMARK(M,C,K,dt,p,u0,v0,a0,[alpha delta]) applies the Newmark
%   method for the calculation of the nodal displacements u, velocities v and
%   accelerations a of the dynamic system with the system matrices M, C and K
%   due to the excitation p.
%
%   M    Mass matrix (nDof * nDof)
%   C    Damping matrix (nDof * nDof)
%   K    Stiffness matrix (nDof * nDof)
%   dt   Time step of the integration scheme (1 * 1).  Should be small enough
%        to ensure the stability and the precision of the integration scheme.
%   p    Excitation (nDof * N).  p(:,k) corresponds to time point t(k).
%   u0   Displacements at time point t(1)-dt (nDof * 1).  Defaults to zero.
%   v0   Velocities at time point t(1)-dt (nDof * 1).  Defaults to zero.
%   a0   Accelerations at time point t(1)-dt (nDof * 1).  Defaults to zero.
%   u    Displacements (nDof * N).  u(:,k) corresponds to time point t(k).
%   t    Time axis (1 * N), defined as t = [0:N-1] * dt.

% Mattias Schevenels
% April 2004

% INITIALISATION
nDof=size(M,1);
N=size(p,2);
t=[0:N-1]*dt;

% DEFAULT INPUT ARGUMENTS
alpha=varargin{end}(1);
delta=varargin{end}(2);
if length(varargin)>1, u0=varargin{1}; else u0=zeros(nDof,1); end;
if length(varargin)>2, v0=varargin{2}; else v0=zeros(nDof,1); end;
if length(varargin)>3, a0=varargin{3}; else a0=zeros(nDof,1); end;

% CHECK INPUT ARGUMENTS
error(argdimchk( ...
  M,     {'nDof','nDof'}, ...
  C,     {'nDof','nDof'}, ...
  K,     {'nDof','nDof'}, ...
  dt,    { 1    , 1    }, ...
  p,     {'nDof','N'   }, ...
  u0,    {'nDof', 1    }, ...
  v0,    {'nDof', 1    }, ...
  a0,    {'nDof', 1    }));

% DIRECT INTEGRATION: NEWMARK METHOD
% (Course text, p. 7.20-7.24)
L=inv(M/alpha/dt^2+delta*C/alpha/dt+K);
u=zeros(nDof,N);
v=zeros(nDof,N);
a=zeros(nDof,N);

for k=1:N
  switch k
    case 1
      R=p(:,k)+M*(u0/alpha/dt^2+v0/alpha/dt+a0*(1/2/alpha-1))+...
               C*(u0*delta/alpha/dt+v0*(delta/alpha-1)+a0*dt/2*(delta/alpha-2));
      u(:,k)=L*R;
      a(:,k)=u(:,k)/alpha/dt^2-u0/alpha/dt^2-v0/alpha/dt-a0*(1/2/alpha-1);
      v(:,k)=v0+dt*(1-delta)*a0+dt*delta*a(:,k);
    otherwise
      R=p(:,k)+M*(u(:,k-1)/alpha/dt^2+v(:,k-1)/alpha/dt+a(:,k-1)*(1/2/alpha-1))+...
               C*(u(:,k-1)*delta/alpha/dt+v(:,k-1)*(delta/alpha-1)+a(:,k-1)*dt/2*(delta/alpha-2));
      u(:,k)=L*R;
      a(:,k)=u(:,k)/alpha/dt^2-u(:,k-1)/alpha/dt^2-v(:,k-1)/alpha/dt-a(:,k-1)*(1/2/alpha-1);
      v(:,k)=v(:,k-1)+dt*(1-delta)*a(:,k-1)+dt*delta*a(:,k);
  end
end
