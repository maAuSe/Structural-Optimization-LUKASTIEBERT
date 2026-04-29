function [KeLCS,MeLCS] = kelcs_shell2(r1,phi,L,h,E,nu)

%KELCS_SHELL2   SHELL2 element stiffness matrix in local coordinate system.
%
%   KeLCS = kelcs_shell2(r1,phi,L,h,E,nu) returns the element stiffness matrix
%   in the local coordinate system for a two-node axisymmetric shell element.
%   The global y-axis is assumed to be the axis of symmetry.
%
%   r1         Radial coordinate of node 1
%   phi        Slope with respect to the xz-plane
%   h          Shell thickness
%   E          Young's modulus
%   nu         Poisson coefficient
%   KeLCS      Element stiffness matrix (6 * 6)

% Mattias Schevenels
% April 2020

% This element is based on Zienkiwicz vol. 2 p. 245.  It is similar to SHELL61
% in ANSYS without extra shape functions, i.e. SHELL61 with KEYOPT(3) = 1.
% However, in ANSYS, a different elasticity matrix D is used.  Reverse
% engineering indicates that it is defined as follows in ANSYS:
%
%   r = interp1([-1,1],[r1;r2],xi);
%   e = -h^2/12*sin(phi)/r;
%   D = E*h/(1-nu^2)*[   1,   nu,         e,      nu*e;
%                       nu,    1,      nu*e,         e;
%                        e, nu*e,    h^2/12, nu*h^2/12;
%                     nu*e,    e, nu*h^2/12,    h^2/12];
%
% where e is the eccentricity of the resultant of constant normal stress in the
% meridional direction.  The first and the third row of the elasticity matrix D
% can be derived by not neglecting the terms z/r_y and z/r_z in the relations
% between member forces and stresses (equation (1-1) in Billington), see also
% the section on the SHELL61 element in the ANSYS Theory Reference for the
% Mechanical APDL and Mechanical Applications.  It is not clear how to derive
% the second and the fourth row.
%
% In Stabil, Zienkiwicz and Billington are followed: the terms z/r_y and z/r_z
% are neglected, leading to e = 0.

% RADIAL COORDINATE OF NODE 2
r2 = r1+L*cos(phi);

% GAUSSIAN QUADRATURE POINTS AND WEIGHTS
Xi = [-1;0;1]*sqrt(0.6);    % ANSYS uses 3 integration points along length
HXi = [5;8;5]/9;
nXi = length(Xi);

% INTEGRATION
KeLCS = zeros(6,6);
for iXi = 1:nXi
  xi = Xi(iXi);
  Hxi = HXi(iXi);
  r = interp1([-1,1],[r1;r2],xi);
  % e = -h^2/12*sin(phi)/r; % use this line for exact correspondence with SHELL61 with KEYOPT(3) = 1.
  e = 0;
  D = E*h/(1-nu^2)*[   1,   nu,         e,      nu*e;
                      nu,    1,      nu*e,         e;
                       e, nu*e,    h^2/12, nu*h^2/12;
                    nu*e,    e, nu*h^2/12,    h^2/12];
  B1 = [                      -1/L,                                      0,                                                 0;
        -(cos(phi)*(xi - 1))/(2*r),  -(sin(phi)*(xi - 1)^2*(xi + 2))/(4*r),         -(sin(phi)*(xi - 1)^2*(xi + 1))/(4*r)*L/2;
                                 0,                            -(6*xi)/L^2,                               -(6*xi - 2)/L^2*L/2;
                                 0, -(2*cos(phi)*((3*xi^2)/4 - 3/4))/(L*r),      (cos(phi)*(- 3*xi^2 + 2*xi + 1))/(2*L*r)*L/2];
  B2 = [                       1/L,                                      0,                                                 0;
         (cos(phi)*(xi + 1))/(2*r),   (sin(phi)*(xi + 1)^2*(xi - 2))/(4*r),         -(sin(phi)*(xi - 1)*(xi + 1)^2)/(4*r)*L/2;
                                 0,                             (6*xi)/L^2,                               -(6*xi + 2)/L^2*L/2;
                                 0,  (2*cos(phi)*((3*xi^2)/4 - 3/4))/(L*r), -(2*cos(phi)*((3*xi^2)/4 + xi/2 - 1/4))/(L*r)*L/2];
  KeLCS(1:3,1:3) = KeLCS(1:3,1:3) + Hxi * pi*r*L*B1'*D*B1;
  KeLCS(1:3,4:6) = KeLCS(1:3,4:6) + Hxi * pi*r*L*B1'*D*B2;
  KeLCS(4:6,1:3) = KeLCS(4:6,1:3) + Hxi * pi*r*L*B2'*D*B1;
  KeLCS(4:6,4:6) = KeLCS(4:6,4:6) + Hxi * pi*r*L*B2'*D*B2;
end

% MAKE ELEMENT STIFFNESS MATRIX PERFECTLY SYMMETRIC
KeLCS = (KeLCS+KeLCS.')/2;

