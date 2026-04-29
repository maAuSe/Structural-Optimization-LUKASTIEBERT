function [K,F,M]=addconstr(Constr,DOF,K,F,varargin)

%ADDCONSTR   Add constraint equations to the stiffness matrix and load vector.
%
%     [K,F]=addconstr(Constr,DOF,K,F)
%   [K,F,M]=addconstr(Constr,DOF,K,[],M)
%   [K,F,M]=addconstr(Constr,DOF,K,F,M)
%   modifies the stiffness matrix, the mass matrix and the load vector according 
%   to the applied constraint equations. The dimensions of the stiffness matrix,
%   the mass matrix and the load vector are kept the same. The resulting 
%   stiffness and mass matrix are not symmetric anymore. This function can be 
%   used as well to apply imposed displacements. 
%
%   Constr     Constraint equation:  
%               Constant=CoefS*SlaveDOF+CoefM1*MasterDOF1+CoefM2*MasterDOF2+...
%              [Constant CoefS SlaveDOF CoefM1 MasterDOF1 CoefM2 MasterDOF2 ...]
%   DOF        Degrees of freedom  (nDOF * 1)
%   K          Stiffness matrix (nDOF * nDOF)
%   F          Load vector  (nDOF * nSteps)
%   M          Mass matrix (nDOF * nDOF)

% David Dooms
% March 2008

% PREPROCESSING
if nargin>4
    M=varargin{1};
end
DOF=DOF(:);
if any(Constr(:,2)==0)
    error('The coefficient of the SlaveDOF is equal to 0 in equation %i.',find(Constr(:,2)==0,1))
end

Const=-Constr(:,1)./Constr(:,2);
SlaveDOF=Constr(:,3);
nEqn=size(Constr,1);

if size(Constr,2)>4
    MasterDOF=Constr(:,5:2:size(Constr,2));
    CoefM=-Constr(:,4:2:(size(Constr,2)-1))./repmat(Constr(:,2),1,size(MasterDOF,2));
end

%% option 1: K not symmetric anymore

nDOF=length(DOF);
remainingDOF=DOF;

for iEqn=1:nEqn
    if isempty(find(abs(remainingDOF-SlaveDOF(iEqn,1))<0.0001))
        error('Degree of freedom %.2f is not a MasterDOF anymore.',SlaveDOF(iEqn,1))
    end
    locSlave=find(abs(DOF-SlaveDOF(iEqn,1))<0.0001);
    Knewrow=zeros(1,nDOF);
    Knewrow(1,locSlave)=1;
    if exist('MasterDOF')
        for idof=1:size(MasterDOF,2)
            if isnan(CoefM(iEqn,idof)) | isinf(CoefM(iEqn,idof))  
                warning('Degree of freedom %.2f is ignored.',MasterDOF(iEqn,idof));
            elseif isempty(find(abs(remainingDOF-MasterDOF(iEqn,idof))<0.0001))
                error('Degree of freedom %.2f is not a MasterDOF anymore.',MasterDOF(iEqn,idof))
            else 
                locdof=find(abs(DOF-MasterDOF(iEqn,idof))<0.0001);
                K(locdof,:)=K(locdof,:)+CoefM(iEqn,idof)*K(locSlave,:);
                K(:,locdof)=K(:,locdof)+CoefM(iEqn,idof)*K(:,locSlave);
                if nargin>4
                    M(locdof,:)=M(locdof,:)+CoefM(iEqn,idof)*M(locSlave,:);
                    M(:,locdof)=M(:,locdof)+CoefM(iEqn,idof)*M(:,locSlave);
                end
                if ~isempty(F)
                    F(locdof,:)=F(locdof,:)+CoefM(iEqn,idof)*F(locSlave,:);
                end
                Knewrow(1,locdof)=-CoefM(iEqn,idof);
            end
        end
    end
    if ~isempty(F)
        F=F+Const(iEqn,1)*repmat(K(:,locSlave),1,size(F,2));
        F(locSlave,:)=-Const(iEqn,1);
    end
    K(:,locSlave)=zeros(nDOF,1);
    K(locSlave,:)=Knewrow;
    if nargin>4
        M(locSlave,:)=zeros(1,nDOF);
        M(:,locSlave)=zeros(nDOF,1);
    end
    remainingDOF=unselectdof(remainingDOF,SlaveDOF(iEqn,1))*remainingDOF;
end
