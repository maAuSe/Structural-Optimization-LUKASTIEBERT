function SeLCS = selcs_plane3(Node,Section,Material,UeLCS,Options)

%SELCS_PLANE3   Compute the element stresses for a plane3 element.
%
%   [SeLCS] = selcs_plane3(Node,Section,Material,UeGCS,Options)
%   computes the element stresses in the
%   local coordinate system for the plane3 element.
%
%   Node       Node definitions    [x y z] (4 * 3)
%   Section    Section definition  [h]  (only used in plane stress)
%   Material   Material definition [E nu rho]
%   UeLCS      Displacements (6 * nSteps)
%   Options    Element options            {Option1 Option2 ...}
%   SeLCS      Element stresses in LCS in corner nodes IJKL
%              (9 * nTimeSteps) [sxx syy sxy]
%
%   See also ELEMSTRESS, SE_PLANE3.

% Stijn François
% 2016

% Constitutive matrix
C=cmat_isotropic(Options.problem,Section,Material);

% Triangle shape function
X=Node(:,1);
Y=Node(:,2);
b1=Y(2)-Y(3); b2=Y(3)-Y(1); b3=Y(1)-Y(2);
c1=X(3)-X(2); c2=X(1)-X(3); c3=X(2)-X(1);

% Element area
Delta=0.5*det([1 X(1) Y(1) 
               1 X(2) Y(2) 
               1 X(3) Y(3)]);

% Shape function derivatives
Be=1/(2*Delta)*[b1  0  b2   0  b3   0 
                0  c1  0  c2   0  c3
                c1  b1 c2  b2  c3  b3];

SeLCS=C*Be*UeLCS;
end

