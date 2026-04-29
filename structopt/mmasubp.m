function [f0,varargout] = mmasubp(subp,varargin)

%MMASUBP   Evaluate the MMA subproblem approximation for specified values x.
%
%   [f0,f1,f2,...] = MMASUBP(subp,x1,x2,...) evaluates the MMA subproblem
%   approximation described by the structure array subp at the values (x1,x2,...).
%
%   INPUT ARGUMENTS
%
%   subp    Structure array containing the parameters describing the MMA
%           subproblem, as returned by the function MMA.
%   x1      Values of the variable x1 (p * q * r * ...).
%   x2      Values of the variable x2 (p * q * r * ...).
%
%   OUTPUT ARGUMENTS
%
%   f0      MMA approximation of the objective function at x (p * q * r * ...).
%   f1      MMA approximation of the 1st constraint at x (p * q * r * ...).
%   f2      MMA approximation of the 2nd constraint at x (p * q * r * ...).
%

% Mattias Schevenels
% March 2017

% RETRIEVE SUBPROBLEM PARAMETERS
p0 = subp.p0;
q0 = subp.q0;
r0 = subp.r0;
p = full(subp.p);
q = full(subp.q);
r = subp.r;
low = subp.low;
upp = subp.upp;
alpha = subp.alpha;
beta = subp.beta;

% REPLACE INF ASYMPTOTES WITH LARGE NUMBERS
low(isinf(low)) = sign(low(isinf(low))).*(beta(isinf(low))-alpha(isinf(low)))*1e4;
upp(isinf(upp)) = sign(upp(isinf(upp))).*(beta(isinf(upp))-alpha(isinf(upp)))*1e4;

% COLLECT INPUT ARGUMENTS
x = varargin;

% DETERMINE PROBLEM DIMENSIONS
m = size(r,1);
n = size(p0,1);
pqr = size(x{1});

% CHECK INPUT ARGUMENT DIMENSIONS
if length(x)~=n, error('%i input arguments x1,x2,... expected.',n); end
for k=1:n
  checkdim(x{k},pqr,'All input argument x1, x2,... must have the same dimensions (p * q * r * ...).');
end

% EVALUATE MMA APPROXIMATION
f0 = repmat(r0,pqr);
f = {};
for j = 1:m
  f{j} = repmat(r(j),pqr);
end
for k = 1:n
  Ux = upp(k)-x{k};
  xL = x{k}-low(k);
  Ux(Ux==0) = realmin;
  xL(xL==0) = realmin;
  f0 = f0 + p0(k)./Ux + q0(k)./xL;
  if p0(k)~=0, f0(Ux<0) = nan; end
  if q0(k)~=0; f0(xL<0) = nan; end
  for j=1:m
    f{j} = f{j} + p(j,k)./Ux + q(j,k)./xL;
    if p(j,k)~=0, f{j}(Ux<0) = realmax; end
    if q(j,k)~=0, f{j}(xL<0) = realmax; end
  end
end

% OUTPUT ARGUMENTS
varargout = f;

% SUBFUNCTION CHECKDIM - CHECK INPUT ARGUMENT DIMENSIONS
function checkdim(x,dim,err)
if length(size(x))~=length(dim) || any(size(x)~=dim), error(err); end

