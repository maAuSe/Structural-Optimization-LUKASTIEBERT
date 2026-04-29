function [L,K,F,M]=lconstr(Constr,DOF,K,F,varargin)

%LCONSTR   Add linear constraint equations to the stiffness matrix and load vector 
%                                                              using Lagrange multipliers.
%
%   [L,K,F]=lconstr(Constr,DOF,K,F)
%   [L,K,F,M]=lconstr(Constr,DOF,K,[],M)
%   [L,K,F,M]=lconstr(Constr,DOF,K,F,M)
%   modifies the stiffness matrix, the mass matrix and the load vector according 
%   to the applied constraint equations. The dimensions of the stiffness matrix,
%   the mass matrix and the load vector increase with the number of constraints.
%   The resulting stiffness and mass matrix are symmetric.  
%
%   Constr     Constraint equation:  
%               Constant=CoefS*SlaveDOF+CoefM1*MasterDOF1+CoefM2*MasterDOF2+...
%              [Constant CoefS SlaveDOF CoefM1 MasterDOF1 CoefM2 MasterDOF2 ...]
%   DOF        Degrees of freedom  (nDOF * 1)
%   L          Selection matrix for displacements: U = L*(K\F);     (nDOF * (nDOF+nConstr))
%   K          Stiffness matrix (nDOF * nDOF)
%   F          Load vector  (nDOF * nSteps)
%   M          Mass matrix (nDOF * nDOF)

% Miche Jansen
% 2012

% PREPROCESSING
if nargin>4
    M=varargin{1};
end
DOF=DOF(:);
nDOF = length(DOF);

nConstr = size(Constr,1);
Constants = Constr(:,1);
sA = Constr(:,2:2:end-1);
CDOF = Constr(:,3:2:end);
nCDOF= size(CDOF,2);

L = spdiags([ones(nDOF,1);zeros(nConstr,1)],0,nDOF,nDOF+nConstr);
nSteps = size(F,2);
F = [F;repmat(Constants,1,nSteps)];

[LIA,jA] = ismember(CDOF,DOF);
iA = repmat((1:nConstr).',nCDOF);
A = sparse(iA(LIA),jA(LIA),sA(LIA),nConstr,nDOF);
K = [K,A.';A,sparse(nConstr,nConstr)];
if nargin>4 
M = [M,sparse(nDOF,nConstr);sparse(nConstr,nDOF),sparse(nConstr,nConstr)];
else
M = [];    
end

