function x=msupt(omega,xi,t,pm,varargin);

%MSUPT   Modal superposition in the time domain.
%
%   x = MSUPT(omega,xi,t,pm,x1,y1,interp) calculates the modal displacements
%   x(t) of a dynamic system with the eigenfrequencies omega and the modal
%   damping ratios xi due to the modal excitation pm(t), given the initial
%   conditions x1 and y1.
%   The solution of the modal differential equations is performed by means of
%   the piecewise exact method.  Interp can be 'foh' (first order hold) or 'zoh'
%   (zero order hold).  In the first case the excitation is assumed to vary
%   linearly within each time step while in the second case it is assumed to be
%   constant: if t(k) <= t < t(k+1) then p(t) = p(k).
%
%   omega   Eigenfrequencies [rad/s] (nMode * 1)
%   xi      Modal damping ratios [ - ] (nMode * 1) or (1 * 1), constant modal
%           damping is assumed in the (1 * 1) case.
%   t       Time points (1 * N) of the sampling of p, x and t.
%   pm      Modal excitation (nMode * N).
%   x1      Modal displacements at time point t(1) (nMode * 1), defaults to zero.
%   y1      Modal velocities at time point t(1) (nMode * 1), defaults to zero.
%   interp  Interpolation scheme: 'foh' or 'zoh', defaults to 'foh'.
%   x       Modal displacements (nMode * N).

% Mattias Schevenels
% March 2004

% INITIALISATION
nMode=length(omega);          % Number of modes to take into account
N=size(pm,2);                 % Number of samples
omegaD=omega.*sqrt(1-xi.^2);  % Eigenfrequencies of the damped system

% DEFAULT INPUT ARGUMENTS
if length(xi)==1, xi=repmat(xi,nMode,1); end; % Constant modal damping
x1=zeros(nMode,1);                            % Zero initial displacements
y1=zeros(nMode,1);                            % Zero initial velocities
order='foh';                                  % First order hold
for iarg=1:length(varargin)
  if isstr(varargin{iarg})
    order=varargin{iarg};
  elseif length(varargin{iarg})==nMode
    if iarg==1, x1=varargin{iarg}; else, y1=varargin{iarg}; end;
  end
end

% CHECK INPUT ARGUMENTS
error(argdimchk( ...
  omega,{'nMode', 1 }, ...
  xi,   {'nMode', 1 }, ...
  t,    { 1     ,'N'}, ...
  pm,   {'nMode','N'}, ...
  x1,   {'nMode', 1 }, ...
  y1,   {'nMode', 1 }));

% INITIALISE x AND y
x=[x1,zeros(nMode,N-1)];
y=[y1,zeros(nMode,N-1)];

% LOOP OVER ALL TIME POINTS
% (Course text, p. 3.22-3.23)
for k=2:N
  dt=t(k)-t(k-1);

  alpha=pm(:,k-1);
  if strcmpi(order,'foh'), beta=(pm(:,k)-pm(:,k-1))/dt; else, beta=0; end;

  a0=alpha./omega.^2-2*xi.*beta./omega.^3;
  a1=beta./omega.^2;

  A=((y(:,k-1)-a1)+xi.*omega.*(x(:,k-1)-a0))./omegaD;
  B=x(:,k-1)-a0;

  x(:,k)=exp(-xi.*omega*dt).*(A.*sin(omegaD*dt)+B.*cos(omegaD*dt))+a0+a1*dt;
  y(:,k)=-xi.*omega.*exp(-xi.*omega*dt).*(A.*sin(omegaD*dt)+B.*cos(omegaD*dt))+ ...
         exp(-xi.*omega*dt).*(A.*omegaD.*cos(omegaD*dt)-B.*omegaD.*sin(omegaD*dt))+a1;

end
