function DOF=removedof(DOF,seldof)

%REMOVEDOF    Remove DOF with Dirichlet boundary conditions equal to zero.
%
%   DOF=removedof(DOF,seldof) removes the specified Dirichlet boundary 
%   conditions from the degrees of freedom vector.
%
%   DOF        Degrees of freedom  (nDOF * 1)
%   seldof  Dirichlet boundary conditions equal to zero      [NodID.dof]
%
%   See also GETDOF.

% David Dooms
% March 2008

% PREPROCESSING
DOF=DOF(:);
seldof=seldof(:);

if ~ isempty(find(seldof==0.00))
    error('The wild card 0.00 is not allowed')
end

nDOF=length(DOF);
ndof=length(seldof);
indj=ones(1,nDOF);

for idof=1:ndof
    if floor(seldof(idof,1))==0       % wild cards 0.0X
        indjdof=find(abs(rem(DOF,1)-rem(seldof(idof,1),1))<0.0001);
    elseif rem(seldof(idof,1),1)==0   % wild cards X.00
        indjdof=find(abs(floor(DOF)-floor(seldof(idof,1)))<0.0001);
    else                                 % standard case
        indjdof=find(abs(DOF-seldof(idof,1))<0.0001);
    end
    if ~isempty(indjdof)
       indj(indjdof)=0;
    end
    
end

DOF=DOF(find(indj));

