function [FdiagrLCS,loc,Extrema,dFdiagrLCSdx] = fdiagrlcs_beam(ftype,Forces,DLoadLCS,L,Points,dForcesdx,dDLoadLCSdx,dLdx)

%FDIAGRLCS_BEAM   Force diagram for a beam element in LCS.
%
%   [FdiagrLCS,loc,Extrema] = FDIAGRLCS_BEAM(ftype,Forces,DLoadLCS,L,Points) 
%   computes the elements forces at the specified points. The extreme
%   values for an element with a single DLoad are analytically determined.
%   The extreme values for an element with multiple DLoads are calculated
%   in the interpolation points only.
%
%   [FdiagrLCS,loc,Extrema,dFdiagrLCSdx] 
%       = FDIAGRLCS_BEAM(ftype,Forces,DLoadLCS,L,Points,dForcesdx,dDLoadLCSdx,dLdx)
%   additionally computes the derivatives of the element force values with
%   respect to the design variables x.
%
%   ftype           'norm'       Normal force (in the local x-direction)
%                   'sheary'     Shear force in the local y-direction
%                   'shearz'     Shear force in the local z-direction
%                   'momx'       Torsional moment (around the local x-direction)
%                   'momy'       Bending moment around the local y-direction
%                   'momz'       Bending moment around the local z-direction
%   Forces          Element forces in LCS (beam convention)	[N; Vy; Vz; T; My; Mz](12 * nLC)
%   DLoadLCS        Distributed loads in LCS [n1localX; n1localY; n1localZ; ...] (6 * nLC)
%   L               Beam length
%   Points          Points in the local coordinate system	(1 * nPoints)
%   dForcesdx       Element forces in LCS derivatives       (SIZE(Forces) * nVar)
%   dDLoadLCSdx     Distributed loads derivatives           (SIZE(DLoadLCS) * nVar)
%   dLdx            Beam length derivatives                 (1 * nVar)
%   FdiagrLCS       Element forces at the points            (1 * nPoints * nLC)
%   loc             Locations of the extreme values         (nValues * nLC)
%   Extrema         Extreme values                          (nValues * nLC)
%                       loc and Extrema are only calculated when nLC = 1 (for plotting). 
%                       If this is not the case their calculation is omitted for effiency.
%   dFdiagrLCSdx    Element force value derivatives         (1 * nPoints * nLC * nVar)
%
%   See also FDIAGRGCS_BEAM.

% David Dooms, Wouter Dillen
% October 2008, December 2016

if nargin<6, dForcesdx = []; end
if nargin<7, dDLoadLCSdx = []; 
else,        dDLoadLCSdx = permute(dDLoadLCSdx,[1 2 4 3]); end
if nargin<8, dLdx = []; end

nVar = 0;
if nargout>3 && (~isempty(dForcesdx) || ~isempty(dDLoadLCSdx) || ~isempty(dLdx))
    nVar = max([size(dForcesdx,3),size(dDLoadLCSdx,4),length(dLdx)]);
end

if nVar==0 || isempty(dForcesdx), dForcesdx = zeros(size(Forces,1),size(Forces,2),nVar); end
if nVar==0 || isempty(dDLoadLCSdx), dDLoadLCSdx = zeros(size(DLoadLCS,1),size(DLoadLCS,2),size(DLoadLCS,3),nVar); end
if nVar==0 || isempty(dLdx), dLdx = zeros(1,nVar); end


Points=Points(:).';
nPoints=length(Points);
nLC=size(DLoadLCS,2);
nDLoadsOnElement = size(DLoadLCS,3);
partialDLoad = -ones(nDLoadsOnElement,1);

FdiagrLCS = zeros(1,nPoints,nLC);
dFdiagrLCSdx = zeros(1,nPoints,nLC,nVar);


%% 1. COMPUTE FORCE POLYNOMIALS 
% (A) Calculate force polynomials for the distributed loads
iFp = repmat(2,[nPoints,nDLoadsOnElement]); % force polynomial column index
A = zeros(4,nLC,3,nDLoadsOnElement);
dA = zeros(4,nLC,3,nDLoadsOnElement,nVar);
v0=zeros(1,nLC);


for i=1:nDLoadsOnElement
    
    % DLoad categorization
    if ~(nnz(DLoadLCS(1:6,:,i))==0)
        if size(DLoadLCS,1)==8 && ~isnan(sum(sum(DLoadLCS(7:8,:,i))))
            partialDLoad(i)=1;
        else 
            partialDLoad(i)=0;
        end
    end

    % COMPLETE DLOAD: evaluate all load cases at once.
    if partialDLoad(i)==0
        
        % Force polynomials for a DLoad
        switch lower(ftype)
            case 'norm'
                A(:,:,2,i) = [v0;  (DLoadLCS(1,:,i)-DLoadLCS(4,:,i))*L/2;  -DLoadLCS(1,:,i)*L;   v0];
            case 'sheary'
                A(:,:,2,i) = [v0;  (DLoadLCS(2,:,i)-DLoadLCS(5,:,i))*L/2;  -DLoadLCS(2,:,i)*L;   v0];
            case 'shearz'
                A(:,:,2,i) = [v0;  (DLoadLCS(3,:,i)-DLoadLCS(6,:,i))*L/2;  -DLoadLCS(3,:,i)*L;   v0];
            case 'momx'
                A(:,:,2,i) = 0;
            case 'momy'
                A(:,:,2,i) = [(DLoadLCS(3,:,i)-DLoadLCS(6,:,i))*L^2/6;  -DLoadLCS(3,:,i)*L^2/2;   v0;  v0];
            case 'momz'
                A(:,:,2,i) = [(DLoadLCS(2,:,i)-DLoadLCS(5,:,i))*L^2/6;  -DLoadLCS(2,:,i)*L^2/2;   v0;  v0];
            otherwise
                error('Unknown element force.')
        end
        % Force polynomials for a DLoad - derivatives
        for n=1:nVar
            dqadx = dDLoadLCSdx(1:3,:,i,n);
            dqbdx = dDLoadLCSdx(4:6,:,i,n);
            switch lower(ftype)
                case 'norm'
                    dA(:,:,2,i,n) = [v0;  (dqadx(1,:)-dqbdx(1,:))*L/2+(DLoadLCS(1,:,i)-DLoadLCS(4,:,i))*dLdx(n)/2;  -dqadx(1,:)*L-DLoadLCS(1,:,i)*dLdx(n);   v0];
                case 'sheary'
                    dA(:,:,2,i,n) = [v0;  (dqadx(2,:)-dqbdx(2,:))*L/2+(DLoadLCS(2,:,i)-DLoadLCS(5,:,i))*dLdx(n)/2;  -dqadx(2,:)*L-DLoadLCS(2,:,i)*dLdx(n);   v0];
                case 'shearz'
                    dA(:,:,2,i,n) = [v0;  (dqadx(3,:)-dqbdx(3,:))*L/2+(DLoadLCS(3,:,i)-DLoadLCS(6,:,i))*dLdx(n)/2;  -dqadx(3,:)*L-DLoadLCS(3,:,i)*dLdx(n);   v0];
                case 'momx'
                    dA(:,:,2,i,n) = 0;
                case 'momy'
                    dA(:,:,2,i,n) = [(dqadx(3,:)-dqbdx(3,:))*L^2/6+(DLoadLCS(3,:,i)-DLoadLCS(6,:,i))*L*dLdx(n)/3;  -dqadx(3,:)*L^2/2-DLoadLCS(3,:,i)*L*dLdx(n);   v0;  v0];
                case 'momz'
                    dA(:,:,2,i,n) = [(dqadx(2,:)-dqbdx(2,:))*L^2/6+(DLoadLCS(2,:,i)-DLoadLCS(5,:,i))*L*dLdx(n)/3;  -dqadx(2,:)*L^2/2-DLoadLCS(2,:,i)*L*dLdx(n);   v0;  v0];
                otherwise
                    error('Unknown element force.')
            end
        end
        

    % PARTIAL DLOAD: evaluate only non-trivial load cases.
    elseif partialDLoad(i)==1
        
        % j = DLoadLCS column index containing the DLoad information (should be only one)
        j = find(sum(abs(DLoadLCS(1:6,:,i)),1));
        if length(j) > 1
            error('Wrong DLoads structure in dimension 3; make sure to use MULTDLOADS when combining multiple load cases.');
        end
    
        % set starting & ending point     
        a = DLoadLCS(7,j,i);
        b = DLoadLCS(8,j,i);

        % FORCE POLYNOMIALS
        % The magnitude of the member force is expressed in function of the local
        % x-coordinate (specified by the interpolation points). The  vectors Ai
        % contain the force polynomials for the (partial) distributed load. For
        % each member force the following 3 situations are considered:
        %       A1: polynomial for all x*L < a
        %       A2: polynomial for all a < x*L < b
        %       A3: polynomial for all x*L > b
        % C is the slope of the distributed load
        switch lower(ftype)
            case 'norm'
                C = (DLoadLCS(4,j,i)-DLoadLCS(1,j,i))/(b-a);
                A(:,j,1,i) = 0;
                A(:,j,2,i) = [0;  -C*L^2/2;  -(DLoadLCS(1,j,i)-a*C)*L;  -a^2*C/2+a*DLoadLCS(1,j,i)];
                A(:,j,3,i) = [0;  0;  0;  -(DLoadLCS(4,j,i)+DLoadLCS(1,j,i))/2*(b-a)];
            case 'sheary'
                C = (DLoadLCS(5,j,i)-DLoadLCS(2,j,i))/(b-a); 
                A(:,j,1,i) = 0;
                A(:,j,2,i) = [0;  -C*L^2/2;  -(DLoadLCS(2,j,i)-a*C)*L;  -a^2*C/2+a*DLoadLCS(2,j,i)];
                A(:,j,3,i) = [0;  0;  0;  -(DLoadLCS(5,j,i)+DLoadLCS(2,j,i))/2*(b-a)];
            case 'shearz'
                C = (DLoadLCS(6,j,i)-DLoadLCS(3,j,i))/(b-a); 
                A(:,j,1,i) = 0;
                A(:,j,2,i) = [0;  -C*L^2/2;  -(DLoadLCS(3,j,i)-a*C)*L;  -a^2*C/2+a*DLoadLCS(3,j,i)];
                A(:,j,3,i) = [0;  0;  0;  -(DLoadLCS(6,j,i)+DLoadLCS(3,j,i))/2*(b-a)]; 
            case 'momx'
                % torsion DLoad on the entire element
                A(:,j,:,i) = 0;
            case 'momy'
                C = (DLoadLCS(6,j,i)-DLoadLCS(3,j,i))/(b-a);
                A(:,j,1,i) = 0;
                A(:,j,2,i) = [-C*L^3/6;    -(DLoadLCS(3,j,i)-C*a)*L^2/2;   (DLoadLCS(3,j,i)*a-C*a^2/2)*L;  -DLoadLCS(3,j,i)*a^2/2+a^3*C/6];
                A(:,j,3,i) = [0;  0;  (-(DLoadLCS(6,j,i)+DLoadLCS(3,j,i))*(b-a)/2)*L;    DLoadLCS(3,j,i)*(b-a)*(b+a)/2+(DLoadLCS(6,j,i)-DLoadLCS(3,j,i))*(2*b^2-a*b-a^2)/6];
            case 'momz'
                C = (DLoadLCS(5,j,i)-DLoadLCS(2,j,i))/(b-a);
                A(:,j,1,i) = 0;
                A(:,j,2,i) = [-C*L^3/6;    -(DLoadLCS(2,j,i)-C*a)*L^2/2;   (DLoadLCS(2,j,i)*a-C*a^2/2)*L;  -DLoadLCS(2,j,i)*a^2/2+a^3*C/6];
                A(:,j,3,i) = [0;  0;  (-(DLoadLCS(5,j,i)+DLoadLCS(2,j,i))*(b-a)/2)*L;    DLoadLCS(2,j,i)*(b-a)*(b+a)/2+(DLoadLCS(5,j,i)-DLoadLCS(2,j,i))*(2*b^2-a*b-a^2)/6];
            otherwise
                error('Unknown element force.')
        end
        for n=1:nVar
            dadx = dDLoadLCSdx(7,j,i,n);
            dbdx = dDLoadLCSdx(8,j,i,n);        
            switch lower(ftype)
                case 'norm'
                    C = (DLoadLCS(4,j,i)-DLoadLCS(1,j,i))/(b-a);
                    dqadx = dDLoadLCSdx(1,j,i,n);
                    dqbdx = dDLoadLCSdx(4,j,i,n);
                    dA(:,j,1,i,n) = 0;
                    dA(:,j,2,i,n) = [0;  -C*L*dLdx(n);  -((dqadx-dadx*C)*L+(DLoadLCS(1,j,i)-C*a)*dLdx(n));  -C*a*dadx+dqadx*a+DLoadLCS(1,j,i)*dadx];
                    dA(:,j,3,i,n) = [0;  0;  0;  -((dqbdx+dqadx)/2*(b-a)+(DLoadLCS(4,j,i)+DLoadLCS(1,j,i))/2*(dbdx-dadx))];                
                case 'sheary'
                    C = (DLoadLCS(5,j,i)-DLoadLCS(2,j,i))/(b-a); 
                    dqadx = dDLoadLCSdx(2,j,i,n);
                    dqbdx = dDLoadLCSdx(5,j,i,n);
                    dA(:,j,1,i,n) = 0;
                    dA(:,j,2,i,n) = [0;  -C*L*dLdx(n);  -((dqadx-dadx*C)*L+(DLoadLCS(2,j,i)-C*a)*dLdx(n));  -C*a*dadx+dqadx*a+DLoadLCS(2,j,i)*dadx];
                    dA(:,j,3,i,n) = [0;  0;  0;  -((dqbdx+dqadx)/2*(b-a)+(DLoadLCS(5,j,i)+DLoadLCS(2,j,i))/2*(dbdx-dadx))];                
                case 'shearz'
                    C = (DLoadLCS(6,j,i)-DLoadLCS(3,j,i))/(b-a); 
                    dqadx = dDLoadLCSdx(3,j,i,n);
                    dqbdx = dDLoadLCSdx(6,j,i,n);
                    dA(:,j,1,i,n) = 0;
                    dA(:,j,2,i,n) = [0;  -C*L*dLdx(n);  -((dqadx-dadx*C)*L+(DLoadLCS(3,j,i)-C*a)*dLdx(n));  -C*a*dadx+dqadx*a+DLoadLCS(3,j,i)*dadx];
                    dA(:,j,3,i,n) = [0;  0;  0;  -((dqbdx+dqadx)/2*(b-a)+(DLoadLCS(6,j,i)+DLoadLCS(3,j,i))/2*(dbdx-dadx))];                
                case 'momx'
                    % torsion DLoad on the entire element
                    dA(:,j,:,i,n) = 0;
                case 'momy'
                    C = (DLoadLCS(6,j,i)-DLoadLCS(3,j,i))/(b-a);
                    dqadx = dDLoadLCSdx(3,j,i,n);
                    dqbdx = dDLoadLCSdx(6,j,i,n);
                    dA(:,j,1,i,n) = 0;
                    dA(:,j,2,i,n) = [-C/6*3*L^2*dLdx(n);  -((dqadx-dadx*C)*L^2/2+(DLoadLCS(3,j,i)-C*a)*L*dLdx(n));  (dqadx*a+DLoadLCS(3,j,i)*dadx-C*a*dadx)*L+(DLoadLCS(3,j,i)*a-C/2*a^2)*dLdx(n);  -a*dadx*DLoadLCS(3,j,i)-a^2/2*dqadx+C/6*3*a^2*dadx];
                    dA(:,j,3,i,n) = [0;  0;  (-(dqbdx+dqadx)*(b-a)/2-(DLoadLCS(6,j,i)+DLoadLCS(3,j,i))*(dbdx-dadx)/2)*L-(DLoadLCS(6,j,i)+DLoadLCS(3,j,i))*(b-a)/2*dLdx(n);  dqadx*(b^2-a^2)/2+DLoadLCS(3,j,i)*(2*b*dbdx-2*a*dadx)/2+(dqbdx-dqadx)*(2*b^2-a*b-a^2)/6+(DLoadLCS(6,j,i)-DLoadLCS(3,j,i))*(4*b*dbdx-dadx*b-a*dbdx-2*a*dadx)/6];
                case 'momz'
                    C = (DLoadLCS(5,j,i)-DLoadLCS(2,j,i))/(b-a);
                    dqadx = dDLoadLCSdx(2,j,i,n);
                    dqbdx = dDLoadLCSdx(5,j,i,n);
                    dA(:,j,1,i,n) = 0;
                    dA(:,j,2,i,n) = [-C/6*3*L^2*dLdx(n);  -((dqadx-dadx*C)*L^2/2+(DLoadLCS(2,j,i)-C*a)*L*dLdx(n));  (dqadx*a+DLoadLCS(2,j,i)*dadx-C*a*dadx)*L+(DLoadLCS(2,j,i)*a-C/2*a^2)*dLdx(n);  -a*dadx*DLoadLCS(2,j,i)-a^2/2*dqadx+C/6*3*a^2*dadx];
                    dA(:,j,3,i,n) = [0;  0;  (-(dqbdx+dqadx)*(b-a)/2-(DLoadLCS(5,j,i)+DLoadLCS(2,j,i))*(dbdx-dadx)/2)*L-(DLoadLCS(5,j,i)+DLoadLCS(2,j,i))*(b-a)/2*dLdx(n);  dqadx*(b^2-a^2)/2+DLoadLCS(5,j,i)*(2*b*dbdx-2*a*dadx)/2+(dqbdx-dqadx)*(2*b^2-a*b-a^2)/6+(DLoadLCS(5,j,i)-DLoadLCS(2,j,i))*(4*b*dbdx-dadx*b-a*dbdx-2*a*dadx)/6];
                otherwise
                    error('Unknown element force.')
            end
        end
        % find force polynomial column indices corresponding to the interpolation points
        iFp(Points<=(a/L),i) = 1; % before
        iFp(Points>=(b/L),i) = 3; % after
    end
end 

% (B) Calculate force polynomials for the member end forces 
switch lower(ftype)
    case 'norm'
        B=[v0;  v0;  v0;   Forces(1,:)];
    case 'sheary'
        B=[v0;  v0;  v0;   Forces(2,:)];
    case 'shearz'
        B=[v0;  v0;  v0;   Forces(3,:)];
    case 'momx'
        B=[v0;  v0;  v0;   Forces(4,:)];
    case 'momy'
        B=[v0;  v0;  Forces(3,:)*L;  Forces(5,:)];
    case 'momz'
        B=[v0;  v0;  Forces(2,:)*L;  Forces(6,:)];
    otherwise
        error('Unknown element force.')
end
dB=zeros(size(B,1),size(B,2),nVar);
for n=1:nVar
    switch lower(ftype)
        case 'norm'
            dB(:,:,n)=[v0;  v0;  v0;   dForcesdx(1,:,n)];
        case 'sheary'
            dB(:,:,n)=[v0;  v0;  v0;   dForcesdx(2,:,n)];
        case 'shearz'
            dB(:,:,n)=[v0;  v0;  v0;   dForcesdx(3,:,n)];
        case 'momx'
            dB(:,:,n)=[v0;  v0;  v0;   dForcesdx(4,:,n)];
        case 'momy'
            dB(:,:,n)=[v0;  v0;  dForcesdx(3,:,n)*L+Forces(3,:)*dLdx(n);  dForcesdx(5,:,n)];
        case 'momz'
            dB(:,:,n)=[v0;  v0;  dForcesdx(2,:,n)*L+Forces(2,:)*dLdx(n);  dForcesdx(6,:,n)];
        otherwise
            error('Unknown element force.')
    end
end



%% 2. PREPARE OUTPUT
% Assuming that starting and ending point of every DLoad are identical for all load cases
loc = [];
Extrema = []; 
    
% FdiagrLCS
dA = permute(dA,[5 1 2 3 4]);
dB = permute(dB,[3 1 2]); 

for i=1:nDLoadsOnElement    
    if partialDLoad(i)==1
        cases = [1,2,3];
    else
        cases = 2;
    end
    for icase = 1:length(cases)
         position = cases(icase);
         ind = iFp(:,i)==position;
         FdiagrLCS(1,ind,:) = FdiagrLCS(1,ind,:) + permute(multpolyval(A(:,:,position,i).',Points(ind)),[3 2 1]);
         dFdiagrLCSdx(1,ind,:,:) = dFdiagrLCSdx(1,ind,:,:) + permute(multpolyval(dA(:,:,:,position,i),Points(ind)),[4 2 3 1]);
    end
end
FdiagrLCS = FdiagrLCS + permute(multpolyval(B.',Points),[3 2 1]);
dFdiagrLCSdx = dFdiagrLCSdx + permute(multpolyval(dB,Points),[4 2 3 1]);
            
% loc, Extrema
if nLC==1
    j = 1;
    
    for i=nDLoadsOnElement   
        
        if partialDLoad(i)==1
            a = DLoadLCS(7,:,i);
            b = DLoadLCS(8,:,i);
        else
            a = zeros(1,nLC);
            b = ones(1,nLC)*L;
        end
    
        if nDLoadsOnElement==1    
            
            % Calculate extreme values analytically 
            if partialDLoad(i)==1
                cases = [1,2,3];
            else
                cases = 2;
            end
            for icase = 1:length(cases)
                position = cases(icase);
                Fp = A(:,j,position,i) + B(:,j);
                loc_tmp = [0; 1];
                if (length(Fp)-find(Fp~=0,1))>1
                    loc_tmp = [loc_tmp; roots(polyder(Fp))];
                end
                loc_tmp = sort(loc_tmp(loc_tmp>=(a/L) & loc_tmp<=(b/L))); 
                Extrema_tmp = polyval(Fp,loc_tmp);
                %
                loc = [loc; loc_tmp(:)];
                Extrema = [Extrema; Extrema_tmp(:)];
            end       
           
        else

            % Calculate extreme values in the interpolation points only
            loc = [0; 1];
            [~,imax] = max(FdiagrLCS(1,:,j)); imax = Points(imax(1));
            [~,imin] = min(FdiagrLCS(1,:,j)); imin = Points(imin(1));
            %
            loc = sort([loc; imax(imax~=0 & imax~=1); imin(imin~=0 & imin~=1)]);
            Extrema = FdiagrLCS(1,ismember(Points,loc),j).';

        end
    end
end

end



function Y = multpolyval(P,X)
    
    % MULTPOLYVAL Evaluate multiple polynomials by matrix multiplication.
    %
    %   Y = MULTPOLYVAL(P,X) returns a matrix where each row contains the
    %   values of the polynomial of the corresponding row in P, evaluated
    %   in X (this is equivalent to row-wise polyval). P is supported up to 
    %   3 array dimensions.
    
    if ~isempty(P)
        
        n = size(P,2);
    
        values = repmat(X(:).',[n,1]);
        for i=1:n
            values(i,:) = values(i,:).^(n-i);
        end
    
        Y = zeros(size(P,1),size(values,2),size(P,3));
        for p=1:size(P,3)
            Y(:,:,p) = P(:,:,p)*values;
        end
        
    else
        
        Y = zeros(size(P,1),length(X(:).'),size(P,3));
        
    end
    
end
