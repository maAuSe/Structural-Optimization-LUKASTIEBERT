%EX3A
%
%   Size and shape optimization of a funicular truss for minimum compliance
%   subject to a volume constraint.

% Mattias Schevenels
% March 2021

% NODES [iNode x y z]
Nodes = [1  0    0    0;
         2  2.5 -0.25 0;
         3  5   -0.5  0;
         4  7.5 -0.75 0;
         5 10   -1    0;
         6 12.5 -0.75 0;
         7 15   -0.5  0;
         8 17.5 -0.25 0;
         9 20    0    0;
        10  2.5  0.5  0;
        11  5    1    0;
        12  7.5  1.5  0;
        13 10    2    0;
        14 12.5  1.5  0;
        15 15    1    0;
        16 17.5  0.5  0;
        17 22.5  7.5  0];

% ELEMENT TYPES {iType type}
Types = {1 'beam'};

% SECTIONS [iSec A ky kz Ixx Iyy Izz]
h = 0.2;
Sections = [1 h^2 inf inf 0 0 h^4/12 h/2 h/2 h/2 h/2];

% MATERIALS [iMat E nu]
Materials = [1 30e9 0.25];

% ELEMENTS [iElt iType iSec iMat iNode1 iNode2 iNode3]
Elements = [1 1 1 1  1  2 17;
            2 1 1 1  2  3 17;
            3 1 1 1  3  4 17;
            4 1 1 1  4  5 17;
            5 1 1 1  5  6 17;
            6 1 1 1  6  7 17;
            7 1 1 1  7  8 17;
            8 1 1 1  8  9 17;
            9 1 1 1  1 10 17;
           10 1 1 1 10 11 17;
           11 1 1 1 11 12 17;
           12 1 1 1 12 13 17;
           13 1 1 1 13 14 17;
           14 1 1 1 14 15 17;
           15 1 1 1 15 16 17;
           16 1 1 1 16  9 17;
           17 1 1 1  2 10 17;
           18 1 1 1  3 11 17;
           19 1 1 1  4 12 17;
           20 1 1 1  5 13 17;
           21 1 1 1  6 14 17;
           22 1 1 1  7 15 17;
           23 1 1 1  8 16 17];

% DEGREES OF FREEDOM
DOF = getdof(Elements,Types);
DOF = removedof(DOF,[0.03,0.04,0.05,1.01,1.02,9.02]);

% MAXIMUM ALLOWABLE VOLUME
Vmax = 2;                         % Maximum volume

% LOAD (PHYSICAL AND ADJOINT)
F = nodalvalues(DOF,[10.02 11.02 12.02 13.02 14.02 15.02 16.02],[-15e4 -15e4 -15e4 -15e4 -15e4 -15e4 -15e4]);

% DESIGN VARIABLES AND BOX CONSTRAINTS
n = 4;                            % Number of design variables
x = [-0.25;-0.5;-0.75;0.2];       % Initial design
xmin = [-2;-2;-2;0.01];
xmax = [0;0;0;1];

% OPTIMIZATION LOOP
iter = 1;                         % Iteration counter
itermax = 500;                    % Maximum number of iterations
change = inf;                     % (Fictitious) initial design change
tol = 1e-4;                       % Convergence threshold
mmaparams = [];                   % Internal MMA parameters
while true

  % UPDATE NODES AND SECTIONS MATRICES
  Nodes([2,3,4],3) = x(1:3)';
  Nodes([8,7,6],3) = x(1:3)';
  Sections(:,2) = x(4)^2;
  Sections(:,7) = x(4)^4/12;
  Sections(:,8:11) = x(4)/2;

  % NODE DEFINITION DERIVATIVES
  dNodesdx = zeros([size(Nodes),n]);
  dNodesdx(2,3,1) = 1;
  dNodesdx(3,3,2) = 1;
  dNodesdx(4,3,3) = 1;
  dNodesdx(6,3,3) = 1;
  dNodesdx(7,3,2) = 1;
  dNodesdx(8,3,1) = 1;

  % SECTION DEFINITION DERIVATIVES
  dSectionsdx = zeros([size(Sections),n]);
  dSectionsdx(:,2,4) = 2*x(4);
  dSectionsdx(:,7,4) = x(4)^3/3;
  dSectionsdx(:,8:11,4) = 0.5;

  % OBJECTIVE FUNCTION AND SENSITIVITIES
  [K,~,dKdx] = asmkm(Nodes,Elements,Types,Sections,Materials,DOF,dNodesdx,dSectionsdx);
  U = K\F;
  C = F'*U;                     % Compliance
  dCdx = zeros(1,n);              % Compliance sensitivity
  for k = 1:n
    dCdx(k) = -U'*dKdx{k}*U;
  end
  f0 = C/1000;
  df0dx = dCdx/1000;

  % CONSTRAINT AND SENSITIVITIES
  [Ve,dVedx] = elemvolumes(Nodes,Elements,Types,Sections,dNodesdx,dSectionsdx);
  V = sum(Ve);
  dVdx = sum(dVedx);
  f = V/Vmax-1;
  dfdx = dVdx/Vmax;

  % GCMMA UPDATE
  [xnew,y,z,lmult,mmaparams,subp,change,history] = gcmma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams);

  % PLOT CURRENT DESIGN
  if ~isinf(change)
    Forces = elemforces(Nodes,Elements,Types,Sections,Materials,DOF,U);
    plotforc('momz',Nodes,Elements,Types,Forces,'Values','off','ForcScal',1e-5); hold('on');
    plotelem(Nodes,Elements,Types,'Numbering','off','GCS','off','LineWidth',2); hold('off');
    axis([0 20 -8 8]);
    drawnow;
  end

  % CHECK CONVERGENCE CRITERIA AND MOVE ON TO NEXT ITERATION
  if change>tol && (iter<itermax || isinf(change))
    iter = iter+1;
    x = xnew;
  else
    break;
  end

end

