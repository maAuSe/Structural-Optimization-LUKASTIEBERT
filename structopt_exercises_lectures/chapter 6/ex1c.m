%EX1C
%
%   Adjoint shape sensitivity analysis of a truss structure.

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

% PHYSICAL AND ADJOINT LOAD VECTORS
Fp = nodalvalues(DOF,[1.02 2.02 3.02 4.02 5.02],[-0.5e4 -1e4 -1e4 -1e4 -1e4]);
Fa = L/umax;
F = [Fp Fa];

% NUMBER OF DESIGN VARIABLES
n = 5;

% NODE DEFINITION DERIVATIVES
dNodesdx = zeros([size(Nodes),n]);
dNodesdx(7,3,1) = 1;
dNodesdx(8,3,2) = 1;
dNodesdx(9,3,3) = 1;
dNodesdx(10,3,4) = 1;
dNodesdx(11,3,5) = 1;

% SECTION DEFINITION DERIVATIVES
dSectionsdx = zeros([size(Sections),n]);

% EVALUATE DISPLACEMENT CONSTRAINT AND ITS SENSITIVITIES
[K,~,dKdx] = asmkm(Nodes,Elements,Types,Sections,Materials,DOF,dNodesdx,dSectionsdx);
U = K\F;
Up = U(:,1);                                             % Physical displacement vector
Ua = U(:,2);                                             % Adjoint displacement vector
f = L'*Up/umax-1;                                        % Displacement constraint
dfdx = zeros(1,n);                                       % Constraint sensitivities
for k = 1:n
  dfdx(k) = -Ua'*dKdx{k}*Up;
end

% DISPLAY RESULTS
fprintf('\n');
fprintf('                 df/dx1          df/dx2          df/dx3          df/dx4         df/dx5\n\n');
fprintf('          %15f %15f %15f %15f %15f\n',dfdx);
fprintf('\n');
