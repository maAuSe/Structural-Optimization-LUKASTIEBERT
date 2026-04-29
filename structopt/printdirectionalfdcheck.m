function printdirectionalfdcheck(h,df0dxp,dfdxp,df0dxp_fd,dfdxp_fd)

%PRINTDIRECTIONALFDCHECK   Print directional finite difference check results.
%
%   PRINTDIRECTIONALFDCHECK(h,df0dxp,dfdxp,df0dxp_fd,dfdxp_fd) prints a table 
%   comparing analytical directional derivatives df_i/dx*p with finite 
%   difference approximations df_i/h for a range of step sizes h.
%
%   For each function f_i, the table shows:
%     - Analytical directional derivative df_i/dx*p
%     - Finite difference approximation (f_i(x+h*p)-f_i(x))/h
%     - Absolute difference |diff|
%
%   INPUT ARGUMENTS
%
%   h          Vector of step sizes (1 * nFDstep).
%   df0dxp     Analytical directional derivative df0/dx*p for the objective (1 * 1).
%   dfdxp      Analytical directional derivatives df/dx*p for the constraints (m * 1).
%   df0dxp_fd  Finite difference directional derivatives (f0(x+h*p)-f0(x))/h for 
%              the objective (1 * nFDstep)
%   dfdxp_fd   Finite difference directional derivatives (f(x+h*p)-f(x))/h for 
%              the constraints (m * nFDstep).

% Mattias Schevenels
% January 2026

h = h(:).';         

labels = [compose('df%idxp',0:size(dfdxp_fd,1)); compose('df%idxp_fd',0:size(dfdxp_fd,1))];

FDdata = [h;reshape(permute(cat(3,repmat([df0dxp;dfdxp],1,numel(h)),[df0dxp_fd;dfdxp_fd],abs(repmat([df0dxp;dfdxp],1,numel(h))-[df0dxp_fd;dfdxp_fd])),[3 1 2]),[],numel(h))];

fprintf('\n    h     ');
fprintf('          %-12s%-12s|diff|',labels{:}); fprintf('\n\n');
fprintf(['    %.4e' repmat('     % .4e % .4e % .4e',1,size(dfdxp_fd,1)+1) '\n'],FDdata);
fprintf('\n');
