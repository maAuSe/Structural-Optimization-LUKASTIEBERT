function [hx,hf0,hf] = plotmmahist(history)

%PLOTMMAHIST   Plot MMA convergence history.
%
%   PLOTMMAHIST(history) plots the convergence history of MMA.
%
%   INPUT PARAMETERS
%
%   history.iter   Iteration counter (1 * k).
%   history.x      Values of the variables x in each iteration (n * k).
%   history.low    Lower asymptotes for x in each iteration (n * k).
%   history.upp    Upper asymptotes for x in each iteration (n * k).
%   history.f0     Objective function in each iteration (1 * k).
%   history.f      Constraints in each iteration (m * k).
%
%   OUTPUT ARGUMENTS
%
%   hx             Handles to axes with convergence plots of the variables x (n * 1).
%   hf0            Handle to axes with convergence plot of the objective function f0 (1 * 1).
%   hf             Handles to axes with convergence plots of the constraints f (m * 1).

% Mattias Schevenels
% March 2017

% RETRIEVE CONVERGENCE HISTORY MATRICES
iter = history.iter;
x = history.x;
if isfield(history,'low'), low = history.low; else low = nan(size(x)); end
if isfield(history,'upp'), upp = history.upp; else upp = nan(size(x)); end
f0 = history.f0;
f = history.f;

% NUMBER OF CONSTRAINTS, DESIGN VARIABLES, ITERATIONS
m = size(f,1);
n = size(x,1);
k = size(x,2);

% TRANSPOSE SINGLE-ROW INPUT ARGUMENTS
if size(iter,2)==1, iter = iter.'; end
if size(f0,2)==1, f0 = f0.'; end

% CHECK INPUT ARGUMENT DIMENSIONS
checkdim(iter,[1 k],'Input argument history.iter must have dimensions (1 * k).')
checkdim(x,[n k],'Input argument history.x must have dimensions (n * k).');
checkdim(low,[n k],'Input argument history.low must have dimensions (n * k) where [n,k] = SIZE(history.x).');
checkdim(upp,[n k],'Input argument history.upp must have dimensions (n * k) where [n,k] = SIZE(history.x).');
checkdim(f0,[1 k],'Input argument history.f0 must have dimensions (1 * k) where k = SIZE(history.x,2).');
checkdim(f,[m k],'Input argument history.f must have dimensions (m * k) where k = SIZE(history.x,2).');

% PLOT EVOLUTION OF X
hx = zeros(n,1);
for j = 1:n
  figure;
  hx(j) = gca;
  colors = get(gca,'ColorOrder');
  blue = colors(1,:);
  red = colors(2,:);
  plot(iter,upp(j,:),'.-','Color',red);
  hold('on');
  plot(iter,x(j,:),'.-','Color',blue);
  plot(iter,low(j,:),'.-','Color',red);
  hold('off');
  xlabel('Iteration');
  title(sprintf('Variable x_%i',j));
  if ~all(isnan(low(:))) && ~all(isnan(upp(:)))
    legend({sprintf('U_%i',j),sprintf('x_%i',j),sprintf('L_%i',j)},'Location','northeastoutside');
  end
end

% PLOT EVOLTUTION OF F0
figure;
hf0 = gca;
plot(iter,f0,'.-');
title('Objective function f_0(x)');
xlabel('Iteration');

% PLOT EVOLUTION OF F
hf = zeros(m,1);
for j = 1:m
  figure;
  hf(j) = gca;
  plot(iter,f(j,:),'.-');
  xlabel('Iteration');
  title(sprintf('Constraint f_%i(x)',j));
end

% DO NOT RETURN UNWANTED OUTPUT ARGUMENTS
if nargout == 0, clear('hx'); end

% SUBFUNCTION CHECKDIM - CHECK INPUT ARGUMENT DIMENSIONS
function checkdim(x,dim,err)
if length(size(x))~=length(dim) || any(size(x)~=dim), error(err); end
