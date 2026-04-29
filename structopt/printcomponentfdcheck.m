function printcomponentfdcheck(h,df0dx,dfdx,df0dx_fd,dfdx_fd,FDvar)

%PRINTCOMPONENTFDCHECK   Print component-wise finite difference check results.
%
%   PRINTCOMPONENTFDCHECK(h,df0dx,dfdx,df0dx_fd,dfdx_fd) prints a table comparing
%   analytical derivatives df_i/dx_j with finite difference approximations
%   for a range of step sizes h.
%
%   PRINTCOMPONENTFDCHECK(...,FDvar) uses the indices FDvar to label the
%   component indices x_j in the table. This argument is only used for
%   labeling; the data arrays are assumed to already correspond to the
%   selected components.
%
%   For each function f_i and each checked component x_j, the table shows:
%     - Analytical derivative df_i/dx_j
%     - Finite difference approximation (f_i(x+h*e_j)-f_i(x))/h
%     - Absolute difference |diff|
%
%   INPUT ARGUMENTS
%
%   h         Vector of step sizes (1 * nFDstep).
%   df0dx     Analytical derivatives df0/dx_j for the objective (1 * nFDvar).
%   dfdx      Analytical derivatives df/dx_j for the constraints (m * nFDvar).
%   df0dx_fd  Finite difference derivatives for the objective (1 * nFDvar * nFDstep).
%   dfdx_fd   Finite difference derivatives for the constraints (m * nFDvar * nFDstep).
%   FDvar     Indices of the checked design variables, used to label the columns
%             (x_j indices) in the printed table (1 * nFDvar). Default: 1:nFDvar.

% Mattias Schevenels
% January 2026

h = h(:).';

m = size(dfdx,1);
nFDvar = size(df0dx_fd,2);

if nargin<6 || isempty(FDvar)
  FDvar = 1:nFDvar;
end

I = repmat((0:m).',1,nFDvar);
J = repmat(reshape(FDvar,1,[]),m+1,1);

ana = "df" + string(I) + "dx" + string(J);
fd  = ana + "_fd";

labels = reshape(permute(cat(3,ana,fd),[3 2 1]),1,[]);

Fana = repmat(cat(1,reshape(df0dx,1,[]),dfdx),1,1,numel(h));   % (m+1) * nFDvar * nFDstep
Ffd  = cat(1,df0dx_fd,dfdx_fd);                                % (m+1) * nFDvar * nFDstep

FDdata = [h;reshape(permute(cat(4,Fana,Ffd,abs(Fana-Ffd)),[4 2 1 3]),[],numel(h))];

fprintf('\n    h     ');
fprintf('          %-12s%-12s|diff|',labels{:}); fprintf('\n\n');
fprintf(['    %.4e' repmat('     % .4e % .4e % .4e',1,(m+1)*nFDvar) '\n'],FDdata);
fprintf('\n');

