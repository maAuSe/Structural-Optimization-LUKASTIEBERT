function [T,Q0,MasterDOF]=tconstr(Constr,DOF)

%TCONSTR     Return matrices to apply constraint equations.
%
%   [T,Q0,MasterDOF]=tconstr(Constr,DOF)
%   returns matrices to apply constraint equations to the stiffness and mass 
%   matrix and the load vector: Kr=T.'*K*T, Mr=T.'*M*T and Fr=T.'*(F-K*Q0). 
%   The original displacement vector is computed using U=T*Ur+Q0. 
%
%   Constr     Constraint equation:  
%               Constant=CoefS*SlaveDOF+CoefM1*MasterDOF1+CoefM2*MasterDOF2+...
%              [Constant CoefS SlaveDOF CoefM1 MasterDOF1 CoefM2 MasterDOF2 ...]
%   DOF        Degrees of freedom  (nDOF * 1)

% David Dooms
% March 2008

% PREPROCESSING
DOF=DOF(:);
if any(Constr(:,2)==0)
    error('The coefficient of the SlaveDOF is equal to 0 in equation %i.',find(Constr(:,2)==0,1))
end

Q=Constr(:,1)./Constr(:,2);
SlaveDOF=Constr(:,3);
nEqn=size(Constr,1);

if ~(length(unique(round(SlaveDOF*100)))==length(SlaveDOF))
        error('Not all slave degrees of freedom are different.')
end

dof=Constr(:,3:2:size(Constr,2));
Coef=-Constr(:,2:2:(size(Constr,2)-1))./repmat(Constr(:,2),1,size(dof,2));

Q0=zeros(size(DOF));

L=unselectdof(DOF,SlaveDOF);
MasterDOF=L*DOF;

T=unselectdof(DOF,SlaveDOF).';
for iEqn=1:nEqn
    locSlave=find(abs(DOF-SlaveDOF(iEqn,1))<0.0001);
    Q0(locSlave,1)=Q(iEqn,1);
    indexCoef=~(isnan(Coef(iEqn,1:end)) | Coef(iEqn,1:end)==0 | isinf(Coef(iEqn,1:end)));
    indexCoef(1,1)=0;
    [L,indexdof]=selectdof(MasterDOF,dof(iEqn,indexCoef)');
    if ~all(isfinite(indexdof))
        error('One of the master degrees of freedom  in equation %i is not correct e.g. the same as a slave degree of freedom.',iEqn)
    elseif ~(length(indexdof)==sum(indexCoef))
        error('Wild cards are not allowed for the degrees of freedom.')
    end
    T(locSlave,indexdof)=Coef(iEqn,indexCoef);
end
