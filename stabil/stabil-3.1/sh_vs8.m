function [Ni,dNi_dxi,dNi_deta,dNi_dzeta] = sh_vs8(xi,eta,zeta)

%SH_VS8     Shape functions for a volume serendipity element with 8 nodes.
%
%   [Ni,dNi_dxi,dNi_deta,dNi_dzeta] = sh_vs8(xi,eta,zeta) 
%   returns the shape functions and its derivatives to the natural coordinates
%   in the point (xi,eta,zeta).
%
%   xi       Natural coordinate                  (scalar)
%   eta      Natural coordinate                  (scalar)
%   zeta     Natural coordinate                  (scalar)
%   Ni       Shape functions in point (xi,eta)   (8 * 1)
%   dN_dxi   Derivative of Ni to xi              (8 * 1)
%   dN_deta  Derivative of Ni to eta             (8 * 1)
%   dN_dzeta Derivative of Ni to zeta            (8 * 1)
%
%   See also KE_SOLID8.

% nodes
nodx = [-1,1,1,-1,-1,1,1,-1].';
nody = [-1,-1,1,1,-1,-1,1,1].';
nodz = [-1,-1,-1,-1,1,1,1,1].';

Nx = (1+nodx*xi);
Ny = (1+nody*eta);
Nz = (1+nodz*zeta);

Ni = Nx.*Ny.*Nz/8;
dNi_dxi = nodx.*Ny.*Nz/8;
dNi_deta = nody.*Nx.*Nz/8;
dNi_dzeta = nodz.*Nx.*Ny/8;

end