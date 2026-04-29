function F=nodalvalues(DOF,seldof,values)

%NODALVALUES   Construct a vector with the values at the selected DOF.
%
%   F=nodalvalues(DOF,seldof,values)
%   constructs a vector with the values at the selected DOF. This function can 
%   be used to obtain a load vector, initial displacements, velocities or 
%   accelerations.
%
%   DOF        Degrees of freedom  (nDOF * 1)
%   seldof     Selected degrees of freedom    [NodID.dof] (nValues * 1)
%   values     Corresponding values           [Value]     (nValues * nSteps)
%   F          Load vector  (nDOF * nSteps)
%
%   See also ELEMLOADS.

% David Dooms
% March 2008

% PREPROCESSING
if ndims(values)==2 && (size(values,1)==1 || size(values,2)==1)
    values=values(:);
end

L=selectdof(DOF,seldof);
if ~(size(L,1)==size(values,1))
    error('Currently wild cards are not allowed')
end
F=L.'*values;
