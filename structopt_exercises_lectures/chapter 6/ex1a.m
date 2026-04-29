%EX1A
%
%   Numerical shape sensitivity analysis of a truss structure.

% Mattias Schevenels
% March 2021

% NODES [iNode x y z]
Nodes = [1 0 0 0;
         2 1 0 0;
         3 2 0 0;
         4 3 0 0;
         5 4 0 0;
         6 5 0 0;
         7 0 1 0;
         8 1 1 0;
         9 2 1 0;
        10 3 1 0;
        11 4 1 0];

% ELEMENT TYPES {iType type}
Types = {1 'truss'};

% SECTIONS [iSec A]
Sections = [1 0.001];

% MATERIALS [iMat E]
Materials = [1 210e9];

% ELEMENTS [iElt iType iSec iMat iNode1 iNode2]
Elements = [1 1 1 1  1  2;
            2 1 1 1  2  3;
            3 1 1 1  3  4;
            4 1 1 1  4  5;
            5 1 1 1  5  6;
            6 1 1 1  7  8;
            7 1 1 1  8  9;
            8 1 1 1  9 10;
            9 1 1 1 10 11;
           10 1 1 1  1  7;
           11 1 1 1  7  2;
           12 1 1 1  2  8;
           13 1 1 1  8  3;
           14 1 1 1  3  9;
           15 1 1 1  9  4;
           16 1 1 1  4 10;
           17 1 1 1 10  5;
           18 1 1 1  5 11;
           19 1 1 1 11  6];

% DEGREES OF FREEDOM
DOF = getdof(Elements,Types);
DOF = removedof(DOF,[0.03,1.01,7.01,6.02]);

% MAXIMUM ALLOWABLE DISPLACEMENT
umax = 0.005;

% OUTPUT SELECTION VECTOR
L = nodalvalues(DOF,1.02,-1);

% LOAD VECTOR
F = nodalvalues(DOF,[1.02 2.02 3.02 4.02 5.02],[-0.5e4 -1e4 -1e4 -1e4 -1e4]);

% NUMBER OF DESIGN VARIABLES
nVar = 5;

% EVALUATE DISPLACEMENT CONSTRAINT AND ITS SENSITIVITIES
h = logspace(0,-20,21)';                                 % Finite difference step sizes
nStep = length(h);                                       % Number of step sizes to test
K = asmkm(Nodes,Elements,Types,Sections,Materials,DOF);  % Stiffness matrix
U = K\F;                                                 % Displacement vector
f = L'*U/umax-1;                                         % Displacement constraint
dfdx = zeros(nStep,nVar);                                % Constraint sensitivities for different step sizes
for iStep = 1:nStep
  for iVar = 1:nVar
    Nodes1 = Nodes;
    Nodes1(iVar+6,3) = Nodes1(iVar+6,3)+h(iStep);
    K1 = asmkm(Nodes1,Elements,Types,Sections,Materials,DOF);
    U1 = K1\F;
    f1 = L'*U1/umax-1;
    dfdx(iStep,iVar) = (f1-f)/h(iStep);
  end
end

% DISPLAY RESULTS
fprintf('\n');
fprintf('    h            df/dx1          df/dx2          df/dx3          df/dx4         df/dx5\n\n');
fprintf('    %.0e %15f %15f %15f %15f %15f\n',[h,dfdx]');
fprintf('\n');
