function [I_cells, KE_quadrature] = mtop_subcell_stiffness(densityPerX, densityPerY, nu)

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
