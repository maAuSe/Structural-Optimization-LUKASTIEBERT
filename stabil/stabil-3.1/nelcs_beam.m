function NeLCS = nelcs_beam(Points,phi_y,phi_z)

%NELCS_BEAM   Shape functions for a beam element.
%
%   NeLCS = nelcs_beam(Points) 
%   NeLCS = nelcs_beam(Points,phi_y,phi_z) 
%   determines the values of the shape functions in the specified points.
%
%   Points     Points in the local coordinate system (1 * nPoints).
%   phi_y      Shear deformation constant in the y-direction (1 * 1).
%   phi_z      Shear deformation constant in the z-direction (1 * 1).
%   NeLCS      Values (nPoints * 12).
%
%   See also DISP_BEAM.

% David Dooms
% March 2008

Points=Points(:).';

if nargin<2 % Shape functions for a Bernouilli beam
    A=[ 0  0 -1  1;
        2 -3  0  1;
        2 -3  0  1;
        0  0 -1  1;
       -1  2 -1  0;
        1 -2  1  0;
        0  0  1  0;
       -2  3  0  0;
       -2  3  0  0;
        0  0  1  0;
       -1  1  0  0;
        1 -1  0  0;];
else % Shape functions for a Timoshenko beam
    A=[ 0               0                        -1                       1;
        2/(1+phi_y)     -3/(1+phi_y)             -phi_y/(1+phi_y)         1;
        2/(1+phi_z)     -3/(1+phi_z)             -phi_z/(1+phi_z)         1;
        0               0                        -1                       1;
        -1/(1+phi_z)    (4+phi_z)/2/(1+phi_z)    -(2+phi_z)/2/(1+phi_z)   0;
        1/(1+phi_y)     -(4+phi_y)/2/(1+phi_y)   (2+phi_y)/2/(1+phi_y)    0;
        0               0                        1                        0;
        -2/(1+phi_y)    3/(1+phi_y)              phi_y/(1+phi_y)          0;
        -2/(1+phi_z)    3/(1+phi_z)              phi_z/(1+phi_z)          0;
        0               0                        1                        0;
        -1/(1+phi_z)    (2-phi_z)/2/(1+phi_z)    phi_z/2/(1+phi_z)        0;
        1/(1+phi_y)     -(2-phi_y)/2/(1+phi_y)   -phi_y/2/(1+phi_y)       0;];
end
NeLCS=zeros(size(Points,2),12);

for k=1:12
    NeLCS(:,k)=polyval(A(k,:),Points);
end