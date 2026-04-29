%EX1C
%
%   Minimum compliance of a ten-bar truss under one point load using MMA.

%   This code contains the solution of an exercise introduced in the web lectures.
%   Try to solve the exercise independently before inspecting or running the code.

% Mattias Schevenels
% February 2021

% NODES [iNode x y z]
Nodes = [1 18.288 9.144 0;
         2 18.288 0     0;
         3  9.144 9.144 0;
         4  9.144 0     0;
         5  0     9.144 0;
         6  0     0     0];

% ELEMENT TYPES {iType type}
Types = {1 'truss'};

% SECTIONS [iSec A]: use a unit section area to obtain K0
Sections = [1 1;
            2 1;
            3 1;
            4 1;
            5 1;
            6 1;
            7 1;
            8 1;
            9 1;
           10 1];

% MATERIALS [iMat E]: E = 1e7 psi = 68.94e9 N/m^2
Materials = [1 68.948e9];

% ELEMENTS [iElt iType iSec iMat iNode1 iNode2]
Elements = [1 1  1 1 5 3;
            2 1  2 1 3 1;
            3 1  3 1 6 4;
            4 1  4 1 4 2;
            5 1  5 1 4 3;
            6 1  6 1 2 1;
            7 1  7 1 5 4;
            8 1  8 1 6 3;
            9 1  9 1 3 2;
           10 1 10 1 4 1];
nElem = size(Elements,1);

% ELEMENT LENGTHS
x1 = Nodes(Elements(:,5),2);
y1 = Nodes(Elements(:,5),3);
x2 = Nodes(Elements(:,6),2);
y2 = Nodes(Elements(:,6),3);
L = sqrt((x2-x1).^2+(y2-y1).^2);

% BOUNDARY CONDITIONS: 2D analysis; simply supported
DOF = getdof(Elements,Types);
DOF = removedof(DOF,0.03);
DOF = removedof(DOF,[5.00; 6.00]);
nDOF = length(DOF);

% LOADS: 2 point loads of 100 kips = 444822 N
F = nodalvalues(DOF,2.02,-444822);

% DETERMINE CONTRIBUTION K0j OF EACH ELEMENT TO THE STIFFNESS MATRIX
K0 = zeros(nDOF,nDOF,nElem);
for iElem = 1:nElem
  K0(:,:,iElem) = asmkm(Nodes,Elements(iElem,:),Types,Sections,Materials,DOF);
end

% FORMULATION OF THE OPTIMIZATION PROBLEM
Vmax = 0.4;                       % Maximum allowed volume
x = repmat(Vmax/sum(L),nElem,1);  % Initial design
xmin = 1e-4;                      % Minimum section area
xmax = 0.1;                       % Maximum section area

% OPTIMIZATION LOOP
iter = 1;                         % Iteration counter
change = inf;                     % (Fictitious) initial design change
tol = 1e-5;                       % Convergence threshold
mmaparams = [];                   % Internal MMA parameters
while change>tol

  % FINITE ELEMENT ANALYSIS
  K = zeros(nDOF,nDOF);
  for iElem = 1:nElem
    K = K + K0(:,:,iElem)*x(iElem);
  end
  U = K\F;

  % OBJECTIVE FUNCTION AND CONSTRAINT
  C = F'*U;                       % Compliance
  V = L'*x;                       % Volume
  f0 = C/1e4;
  f = V/Vmax-1;

  % SENSITIVITIES
  dCdx = zeros(nElem,1);
  for iElem = 1:nElem
    dCdx(iElem) = -U'*K0(:,:,iElem)*U;
  end
  dVdx = L;
  df0dx = dCdx/1e4;
  dfdx = dVdx/Vmax;

  % MMA UPDATE
  [xnew,y,z,lmult,mmaparams,subp,change,history] = mma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams);

  % PLOT CURRENT DESIGN
  SectionsX = Sections;
  SectionsX(:,2) = x;
  plotelemsec(Nodes,Elements,Types,SectionsX,'Numbering','off','GCS','off');
  drawnow;

  % CHECK CONVERGENCE CRITERION AND MOVE ON TO NEXT ITERATION
  if change>tol
    iter = iter+1;
    x = xnew;
  end

end
