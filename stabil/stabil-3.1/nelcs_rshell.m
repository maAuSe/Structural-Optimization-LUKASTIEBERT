function [Nxe,Nye,Nze]=nelcs_rshell(xi,eta);

% Shape functions in terms of normalized coordinates.  The shape functions
% relating UZ to RX and RY must be multiplied with the half-width of the element
% in the direction Y and X, respectively.

% Mattias Schevenels
% April 2008

xi=xi(:);
eta=eta(:);
N=length(xi);

if length(xi)~=length(eta)
  error('Different number of xi-coordinates and eta-coordinates.');
end

Nxe=1/4*[(1-xi).*(1-eta),zeros(N,5),(1+xi).*(1-eta),zeros(N,5), ...
             (1+xi).*(1+eta),zeros(N,5),(1-xi).*(1+eta),zeros(N,5)];

Nye=1/4*[zeros(N,1),(1-xi).*(1-eta),zeros(N,5),(1+xi).*(1-eta),zeros(N,5), ...
             (1+xi).*(1+eta),zeros(N,5),(1-xi).*(1+eta),zeros(N,4)];

Nze=zeros(N,24);
xii=[-1 1 1 -1];
etai=[-1 -1 1 1];
for n=1:4
  xi0=xi*xii(n);
  eta0=eta*etai(n);
  Nze(:,6*n-3)=1/8*(1+xi0).*(1+eta0).*(2+xi0+eta0-xi.^2-eta.^2);
  Nze(:,6*n-2)=-1/8*(1+xi0).*(1+eta0)*etai(n).*(1-eta.^2);
  Nze(:,6*n-1)=1/8*(1+xi0).*(1+eta0)*xii(n).*(1-xi.^2);
end

