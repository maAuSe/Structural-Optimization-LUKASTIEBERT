function msg=argdimchk(varargin)

%ARGDIMCHK   Validate input argument dimensions.
%
%   msg = ARGDIMCHK(arg1,size1,arg2,size2,...) returns an appropriate error
%   message if the dimensions of arg1,arg2,... do not comply with size1,size2,...
%   respectively.  If they do comply, an empty matrix is returned.
%   size1,size2,... are cell arrays consisting of at least 2 cells.  Each cell
%   corresponds to a dimension of arg1,arg2,... and contains a number (to
%   constrain the dimension explicitly) or a string (to constrain the dimension
%   implicitly).  E.g. the expression:
%
%   ERROR(ARGDIMCHK( ...
%     omega,{'nMode',   1      }, ...
%     Phi,  {'nDof',   'nMode' }, ...
%     xi,   {'nMode',   1      }, ...
%     b,    {'nDof',    1      }, ...
%     q,    { 1,       'nOmega'}, ...
%     Omega,{ 1,       'nOmega'}, ...
%     c,    {'nSelDof','nDof'  }));
%
%   checks if omega,Phi,xi,b,q,Omega,c are all 2-dimensional variables and if
%   SIZE(omega,2)==1           SIZE(Phi,2)==SIZE(omega,1)
%   SIZE(xi,1)==SIZE(omega,1)  SIZE(xi,2)==1
%   SIZE(b,1)==SIZE(Phi,1)     SIZE(b,2)==1
%   SIZE(q,1)==1               SIZE(Omega,1)==1
%   SIZE(Omega,2)==SIZE(q,2)   SIZE(c,2)==SIZE(Phi,1)
%   If not, an appropriate error message is shown.

% Mattias Schevenels
% February 2004

% INITIALISATION
errorLevel=0;    % Defines error type priorities:
                 % - a variable has a wrong number of dimensions      -> 3
                 % - error w.r.t. an explicitly constrained dimension -> 2
                 % - error w.r.t. an implicitly constrained dimension -> 1
                 % - no errors                                        -> 0

msg=[];          % Error message, remains an empty matrix if no errors occur.

symDimValues=[]; % Structure array the fields of which correspond to the
                 % symbolic dimensions, i.e. the names of the implicitly
                 % constrained dimensions.  Each field holds the value of the
                 % first encountered actual dimension linked to the field's
                 % symbolic dimension.

symDimLabels=[]; % Analogous structure array each fields of which holds a label
                 % referring to the first encountered actual dimension linked
                 % to the field's symbolic dimension.

% LOOP OVER ALL ARGUMENTS
for iArg=1:2:nargin
  argName=inputname(iArg);      % Name of the argument in the caller function
  actDims=size(varargin{iArg}); % Actual dimensions of the argument
  actNDim=length(actDims);      % Actual number of dimensions
  expDims=varargin{iArg+1};     % Expected dimensions
  expNDim=length(expDims);      % Expected number of dimensions

  % CHECK NUMBER OF DIMENSIONS
  if (actNDim~=expNDim) & (errorLevel<3)
    msg=['Input argument "' argName '": ' num2str(expNDim) '-dimensional variable expected.'];
    errorLevel=3;
  end

  % CHECK DIMENSION CONSTRAINTS
  for iDim=1:expNDim
    dimLabel=['SIZE(' argName ',' num2str(iDim) ')'];  % Label referring to the current dimension
    actDim=actDims(iDim);                              % Actual value of the current dimension
    expDim=expDims{iDim};                              % Expected value of the current dimension

    % CHECK EXPLICIT DIMENSION CONSTRAINTS
    if ~isstr(expDim)
      if (actDim~=expDim) & (errorLevel<2)
        if (expNDim==2) & (expDim==1) & (iDim==1),
          msg=['Input argument "' argName '": row matrix expected.'];
        elseif (expNDim==2) & (expDim==1) & (iDim==2),
          msg=['Input argument "' argName '": column matrix expected.'];
        else
          msg=['Input argument error: ' dimLabel ' should be ' num2str(expDim) '.'];
        end;
        errorLevel=2;
      end

    % CHECK IMPLICIT DIMENSION CONSTRAINTS
    else
      if isfield(symDimValues,expDim)
        % Check consistency with previously defined symbolically constrained dimensions.
        if (actDim~=getfield(symDimValues,expDim)) & (errorLevel<1)
          msg=['Inconsistent input arguments: ' getfield(symDimLabels,expDim) ' does not match ' dimLabel '.'];
          errorLevel=1;
        end
      else
        % Add the value and the label of the symbolically constrained dimension to the structure arrays.
        symDimValues=setfield(symDimValues,expDim,actDim);
        symDimLabels=setfield(symDimLabels,expDim,dimLabel);
      end
    end

  end

end