function [FdiagrLCS,loc,Extrema] = fdiagrlcs_shell2(ftype,Forces,DLoadLCS,r1,phi,L,h,E,nu,Points)

%FDIAGRLCS_SHELL2   Force diagram for a SHELL2 element in LCS.
%
%   [FdiagrLCS,loc,Extrema] = fdiagrlcs_shell2(ftype,Forces,DLoadLCS,L,Points)
%   computes the elements forces at the specified points. The extreme values are
%   obtained by enumeration.
%
%   ftype      'Nphi'       Normal force (per unit length) in meridional direction
%              'Qphi'       Transverse force (per unit length) in meridional direction
%              'Mphi'       Bending moment (per unit length) in meridional direction
%              'Ntheta'     Normal force (per unit length) in circumferential direction
%              'Mtheta'     Bending moment (per unit length) in circumferential direction
%   Forces     Element forces in LCS (beam convention) [N; Vy; 0; 0; 0; Mz] (12 * 1)
%   DLoadLCS   Distributed loads in LCS [n1localX; n1localY; n1localZ; ...]
%                                                                (6 * 1) or (12 * 1)
%   Points     Points in the local coordinate system (1 * nPoints)
%   FdiagrLCS  Element forces at the points (1 * nPoints)
%   loc        Locations of the extreme values (nValues * 1)
%   Extrema    Extreme values (nValues * 1)

% Mattias Schevenels
% April 2020

% SWITCH FROM BEAM TO ALGEBRAIC CONVENTION
Forces(1:5)=-Forces(1:5);
Forces(12)=-Forces(12);

% PREPROCESSING
r2 = r1+L*cos(phi);
Points = Points(:).';
nPoint = length(Points);
if Points(1)~=0 || Points(end)~=1 || ~issorted(Points)
  error('The Points vector must include 0 and 1 and must be sorted in ascending order.');
end
r = r1+Points*L*cos(phi);

% SUBDIVIDE ELEMENT; ASSEMBLE SUBSYSTEM STIFFNESS MATRIX AND LOAD VECTOR
nElem = nPoint-1;
iK = zeros(6,6,nElem);
jK = zeros(6,6,nElem);
Ke = zeros(6,6,nElem);
[ik,jk] = ndgrid(1:6,1:6);
dload = interp1([0;1],reshape(DLoadLCS,[],2)',Points')';
F0e = zeros(6,nElem);
for iElem = 1:nElem
  r11 = r(iElem);
  L1 = (Points(iElem+1)-Points(iElem))*L;
  Ke(:,:,iElem) = kelcs_shell2(r11,phi,L1,h,E,nu);
  iK(:,:,iElem) = 3*(iElem-1)+ik;
  jK(:,:,iElem) = 3*(iElem-1)+jk;
  F0e(:,iElem) = loadslcs_shell2([dload(:,iElem);dload(:,iElem+1)],L1);
end
F0 = [F0e(1:3,:),[0;0;0]]+[[0;0;0],F0e(4:6,:)]; F0 = F0(:);
K = sparse(iK(:),jK(:),Ke(:));
F1 = zeros(3*nElem+3,1);
F1([1:3,end-2:end]) = Forces([1,2,6,7,8,12]);
F = F0+F1;

% DETERMINE SUBSYSTEM DISPLACEMENTS (RELATIVE TO NODE 2)
U = zeros(3*nElem+3,1);
k = 1:length(U);
if abs(cos(phi))>sqrt(2)/2
  k = k([1:end-2,end]);
else
  k = k([1:end-3,end-1:end]);
end
U(k) = K(k,k)\F(k);
Ue = zeros(6,nElem);
Ue(1:3,:) = reshape(U(1:end-3),3,[]);
Ue(4:6,:) = reshape(U(4:end),3,[]);

% DETERMINE SUBSYSTEM ELEMENT FORCES
Fe = zeros(6,nElem);
for iElem = 1:nElem
  Fe(:,iElem) = Ke(:,:,iElem)*Ue(:,iElem)-F0e(:,iElem);
end

% DETERMINE IN-PLANE MEMBER FORCES
Nphi = [-Fe(1,:) Fe(4,end)]/2/pi./r;
Qphi = [-Fe(2,:) Fe(5,end)]/2/pi./r;
Mphi = [Fe(3,:) -Fe(6,end)]/2/pi./r;

% DETERMINE OUT-OF-PLANE MEMBER FORCES
u = U(1:3:end)';
v = U(2:3:end)';
dvdx = U(3:3:end)';
epsilon1 = zeros(2,nPoint);
epsilon2 = [(u*cos(phi)-v*sin(phi))./r; -dvdx*cos(phi)./r];
sigma1 = [Nphi; Mphi];
sigma2 = zeros(2,nPoint);
for iPoint = 1:nPoint
  % e = -h^2/12*sin(phi)/r(iPoint); % Use this line (and the corresponding line in kelcs_shell2) for compliance with ANSYS.
                                    % However, this leads to counterintuitive Mtheta values, e.g. for a cylinder under internal pressure, Mtheta causes additional tensile stresses in the outer fibers.
                                    % Comparison with ANSYS is not possible, as ANSYS can only compute Mtheta via the stresses.
  e = 0;
  D = E*h/(1-nu^2)*[   1,   nu,         e,      nu*e;
                      nu,    1,      nu*e,         e;
                       e, nu*e,    h^2/12, nu*h^2/12;
                    nu*e,    e, nu*h^2/12,    h^2/12];
  D11 = D([1,3],[1,3]);
  D12 = D([1,3],[2,4]);
  D21 = D([2,4],[1,3]);
  D22 = D([2,4],[2,4]);
  D11(isinf(D11)) = 0;
  epsilon1(:,iPoint) = D11\(sigma1(:,iPoint)-D12*epsilon2(:,iPoint));
  sigma2(:,iPoint) = D21*epsilon1(:,iPoint)+D22*epsilon2(:,iPoint);
end
Ntheta = sigma2(1,:);
Mtheta = sigma2(2,:);

% STRESSES
sNphi = Nphi/h;
sMphiT = 6*Mphi/h^2;
sMphiB = -6*Mphi/h^2;
sNtheta = Ntheta/h;
sMthetaT = 6*Mtheta/h^2;
sMthetaB = -6*Mtheta/h^2;

% ASSIGN OUTPUT DIAGRAM
switch lower(ftype)
  case 'nphi', FdiagrLCS = Nphi;
  case 'qphi', FdiagrLCS = Qphi;
  case 'mphi', FdiagrLCS = Mphi;
  case 'ntheta', FdiagrLCS = Ntheta;
  case 'mtheta', FdiagrLCS = Mtheta;
  case 'snphi', FdiagrLCS = sNphi;
  case 'smphit', FdiagrLCS = sMphiT;
  case 'smphib', FdiagrLCS = sMphiB;
  case 'sntheta', FdiagrLCS = sNtheta;
  case 'smthetat', FdiagrLCS = sMthetaT;
  case 'smthetab', FdiagrLCS = sMthetaB;
  otherwise
    error('Unknown element force.')
end
FdiagrLCS(r==0) = nan;

% DETERMINE EXTREMA
[Min,iMin] = min(FdiagrLCS);
[Max,iMax] = max(FdiagrLCS);
iMin = iMin(~isnan(Min));
iMax = iMax(~isnan(Max));
iExt = unique([1; iMin; iMax; nPoint]);
iExt = iExt(~isnan(FdiagrLCS(iExt)));
loc = Points(iExt)';
Extrema = FdiagrLCS(iExt)';
