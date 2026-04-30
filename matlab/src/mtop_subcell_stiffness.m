function [I_cells, KE_quadrature] = mtop_subcell_stiffness(densityPerX, densityPerY, nu)
%MTOP_SUBCELL_STIFFNESS Per-density-cell stiffness templates for Q4/nN MTOP.
%
% I_cells is an 8 x 8 x N array, with N = densityPerY * densityPerX, such
% that for an analysis element whose density cells carry the values rho_k
% (stored in column-major order on the densityPerY x densityPerX sub-grid),
% the element stiffness matrix used by the multiresolution scheme of
% Nguyen et al. (2010), Eq. (11), is
%
%   K_e = sum_k (E_min + (E_0 - E_min) * rho_k^p) * I_cells(:, :, k),
%
% with the unit-modulus normalization E_0 = 1 used in the lecture code.
% The matrix I_cells(:, :, k) is the contribution of density cell k to the
% midpoint-rule integration of B^T D^0 B over a unit Q4 element. Its
% sum over k therefore approximates the analytical element stiffness
% matrix KE used by the classical baseline; for densityPerX = densityPerY
% = 5 the residual quadrature error is O(h^2) per element, which is the
% expected accuracy of the n25 multiresolution scheme.
%
% The element-local node ordering matches the convention in the Chapter 9
% lecture code and in the 88-line code of Andreassen et al. (2011): node 1
% is the bottom-left corner and the four nodes are listed counter-clockwise.
% Density cells in xPhys follow the standard image array convention
% (sub_y = 1 at the top of the image, sub_y = densityPerY at the bottom),
% so the local y_local of cell (sub_y, sub_x) is mapped to the physical
% top of the element (y_local close to 1) for the smallest sub_y.

C0 = 1 / (1 - nu^2) * [1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];

n_cells = densityPerY * densityPerX;
I_cells = zeros(8, 8, n_cells);
A_phys = 1 / n_cells;

idx = 0;
for sub_x = 1:densityPerX
  for sub_y = 1:densityPerY
    idx = idx + 1;
    x_local = (sub_x - 0.5) / densityPerX;
    y_local = 1 - (sub_y - 0.5) / densityPerY;
    xi = 2 * x_local - 1;
    eta = 2 * y_local - 1;
    B = strain_displacement(xi, eta);
    I_cells(:, :, idx) = B' * C0 * B * A_phys;
  end
end

KE_quadrature = sum(I_cells, 3);

end


function B = strain_displacement(xi, eta)
%STRAIN_DISPLACEMENT Bilinear Q4 strain-displacement matrix on a unit square.
%
% Parent coordinates (xi, eta) lie in [-1, 1]^2 and are mapped onto the
% physical unit square [0, 1]^2 using x_local = (xi + 1) / 2 and
% y_local = (eta + 1) / 2. The four element nodes are ordered as
%   node 1: (0, 0) -> (xi, eta) = (-1, -1)
%   node 2: (1, 0) -> ( 1, -1)
%   node 3: (1, 1) -> ( 1,  1)
%   node 4: (0, 1) -> (-1,  1).

dN_dxi  = [-(1 - eta),  (1 - eta), (1 + eta), -(1 + eta)] / 4;
dN_deta = [-(1 - xi),  -(1 + xi),  (1 + xi),   (1 - xi)] / 4;
dN_dx = 2 * dN_dxi;
dN_dy = 2 * dN_deta;

B = zeros(3, 8);
B(1, 1:2:7) = dN_dx;
B(2, 2:2:8) = dN_dy;
B(3, 1:2:7) = dN_dy;
B(3, 2:2:8) = dN_dx;

end
