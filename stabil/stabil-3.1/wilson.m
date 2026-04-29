function [u,v,a,t]=wilson(M,C,K,dt,p,varargin)

%WILSON   Direct time integration for dynamic systems - Wilson-theta method
%   [u,v,a,t] = WILSON(M,C,K,dt,p,u1,v1,[alpha delta theta]) applies the
%   Wilson-theta method for the calculation of the nodal displacements u,
%   velocities v and accelerations a of the dynamic system with the system
%   matrices M, C and K due to the excitation p.
%
%   M    Mass matrix (nDof * nDof)
%   C    Damping matrix (nDof * nDof)
%   K    Stiffness matrix (nDof * nDof)
%   dt   Time step of the integration scheme (1 * 1).  Should be small enough
%        to ensure the stability and the precision of the integration scheme.
%   p    Excitation (nDof * N).  p(:,k) corresponds to time point t(k).
%   u1   Displacements at time point t(1) (nDof * 1).  Defaults to zero.
%   v1   Velocities at time point t(1) (nDof * 1).  Defaults to zero.
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
theta=varargin{end}(3);
if length(varargin)>1, u1=varargin{1}; else u1=zeros(nDof,1); end;
if length(varargin)>2, v1=varargin{2}; else v1=zeros(nDof,1); end;

% CHECK INPUT ARGUMENTS
error(argdimchk( ...
  M,     {'nDof','nDof'}, ...
  C,     {'nDof','nDof'}, ...
  K,     {'nDof','nDof'}, ...
  dt,    { 1    , 1    }, ...
  p,     {'nDof','N'   }, ...
  u1,    {'nDof', 1    }, ...
  v1,    {'nDof', 1    }));

% DIRECT INTEGRATION: WILSON-THETA METHOD
% (Course text, p. 7.24-7.25)
a1=M\(p(:,1)-K*u1-C*v1);
L=inv(M/alpha/(theta*dt)^2+delta*C/alpha/(theta*dt)+K);
u=zeros(nDof,N);
v=zeros(nDof,N);
a=zeros(nDof,N);
u(:,1)=u1;
v(:,1)=v1;
a(:,1)=a1;

for k=2:N
  pTheta=(1-theta)*p(:,k-1)+theta*p(:,k);
  R=pTheta+ ...
    M*(u(:,k-1)/alpha/(theta*dt)^2+v(:,k-1)/alpha/(theta*dt)+a(:,k-1)*(1/2/alpha-1))+ ...
    C*(u(:,k-1)*delta/alpha/(theta*dt)+v(:,k-1)*(delta/alpha-1)+a(:,k-1)*(theta*dt)/2*(delta/alpha-2));
  uTheta=L*R;
  aTheta=uTheta/alpha/(theta*dt)^2-u(:,k-1)/alpha/(theta*dt)^2-v(:,k-1)/alpha/(theta*dt)-a(:,k-1)*(1/2/alpha-1);
  a(:,k)=a(:,k-1)*(1-1/theta)+aTheta/theta;
  v(:,k)=v(:,k-1)+dt*(1-delta)*a(:,k-1)+dt*delta*a(:,k);
  u(:,k)=u(:,k-1)+dt*v(:,k-1)+dt^2*(1/2-alpha)*a(:,k-1)+dt^2*alpha*a(:,k);
end
