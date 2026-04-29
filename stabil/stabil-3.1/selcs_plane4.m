function SeLCS = selcs_plane4(Node,Section,Material,UeLCS,Options)

%SELCS_PLANE4   Compute the element stresses for a plane4 element.
%
%   [SeLCS] = selcs_plane4(Node,Section,Material,UeGCS,Options)
%   computes the element stresses in the
%   local coordinate system for the plane4 element.
%
%   Node       Node definitions    [x y z] (4 * 3)
%   Section    Section definition  [h]  (only used in plane stress)
%   Material   Material definition [E nu rho]
%   UeLCS      Displacements (8 * nSteps)
%   Options    Element options            {Option1 Option2 ...}
%   SeLCS      Element stresses in LCS in corner nodes IJKL
%              (12 * nTimeSteps) [sxx syy sxy] for plane stress / plane strain problems
%              (16 * nTimeSteps) [sxx stheta syy sxy] for axisymmetric problems
%
%   See also ELEMSTRESS, SE_PLANE4.

% Stijn François
% 2016

% Number of stress components
if strcmpi(Options.problem,'axisym')     % MS 22/04/2020
  n=4;
else
  n=3;
end

% Gaussian points
xg=[-0.577350269189626 -0.577350269189626;
     0.577350269189626 -0.577350269189626;
     0.577350269189626  0.577350269189626;
    -0.577350269189626  0.577350269189626];
Hg=[1 1 1 1];
nGauss=numel(Hg);

if Options.bendingmodes
  K21 = zeros(4,8);
  K22  = zeros(4,4);
  J0=[0.25*[-1 1 1 -1]*Node(:,[1,2]);0.25*[-1 -1 1 1]*Node(:,[1,2])];
  detJ0=det(J0);
end

% Constitutive matrix
C=cmat_isotropic(Options.problem,Section,Material);
CB=zeros(n,8,nGauss);   % MS 22/04/2020: 3 -> n

if Options.bendingmodes, DBex=zeros(3,4,nGauss); end % contribution of bending modes

for iGauss = 1:nGauss
  xi= xg(iGauss,1);
  eta = xg(iGauss,2);

  % Shape functions, Jacobian, ...
  [Ni,dN_dxi,dN_deta] = sh_qs4(xi,eta);

  J = [dN_dxi.' *Node(:,1)  dN_dxi.' *Node(:,2);
       dN_deta.'*Node(:,1)  dN_deta.'*Node(:,2)];
  detJ = det(J);
  dNi = J\[dN_dxi dN_deta].';

  % Shape function derivatives
  if strcmpi(Options.problem,'axisym')     % MS 22/04/2020
    r = Ni.'*Node(:,1);
    B = zeros(4,8);
    B(1,1:2:7) = dNi(1,:);   %epsilon_r
    B(2,1:2:7) = Ni./r;      %epsilon_theta
    B(3,2:2:8) = dNi(2,:);   %epsilon_z
    B(4,1:2:7) = dNi(2,:);   %gamma_rz
    B(4,2:2:8) = dNi(1,:);   %gamma_rz
  else
    B = zeros(3,8);
    B(1,1:2:7) = dNi(1,:);
    B(2,2:2:8) = dNi(2,:);
    B(3,1:2:7) = dNi(2,:);
    B(3,2:2:8) = dNi(1,:);
  end

  CB(:,:,iGauss)=C*B;

  % Bending mode contribution
  if Options.bendingmodes
    Nex = [1-xi^2 -2*xi 0;1-eta^2 0 -2*eta]; % [N5 N6]
    dNex = J0\((detJ0/detJ)*Nex(:,2:3).');
    Bex = zeros(3,4);
    Bex(1,[1 3])=dNex(1,:);
    Bex(2,[2 4])=dNex(2,:);
    Bex(3,[1 3])=dNex(2,:);
    Bex(3,[2 4])=dNex(1,:);
    K21 = K21 + Bex.'*C*B  *Hg(1,iGauss)*detJ;
    K22 = K22 + Bex.'*C*Bex*Hg(1,iGauss)*detJ;
    DBex(:,:,iGauss)=DBex(:,:,iGauss) + C*Bex*Hg(1,iGauss)*detJ;
  end
end


if Options.bendingmodes
  Uex=-K22\K21*UeLCS;  % Internal dofs.
end

SmeLCSg=zeros(4*n,1);                                                     % MS 22/04/2020: 12 -> 4*n
for iGauss = 1:nGauss
  SmeLCSg((1:n)+n*(iGauss-1),:)= CB(:,:,iGauss)*UeLCS;                    % MS 22/04/2020: 3 -> n
  if Options.bendingmodes
    SmeLCSg((1:3)+3*(iGauss-1),:)=  SmeLCSg((1:3)+3*(iGauss-1),:)+DBex(:,:,iGauss)*Uex;
  end
end

g=1/abs(xg(1,1))/2;
extrap = [(1+g)*eye(n) -0.5*eye(n) (1-g)*eye(n) -0.5*eye(n);              % MS 22/04/2020: 3 -> n
          -0.5*eye(n) (1+g)*eye(n) -0.5*eye(n) (1-g)*eye(n);
          (1-g)*eye(n) -0.5*eye(n) (1+g)*eye(n) -0.5*eye(n);
          -0.5*eye(n) (1-g)*eye(n) -0.5*eye(n) (1+g)*eye(n)];
SeLCS = extrap*SmeLCSg;
end
