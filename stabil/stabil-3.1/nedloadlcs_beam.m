function [NeDLoad,dNeDLoaddx] = nedloadlcs_beam(Points,phi_y,phi_z,a,b,L,dadx,dbdx,dLdx)

%NEDLOADLCS_BEAM    Interpolation functions for a distributed load on a beam element.
%
%   NeLCS = NEDLOADLCS_BEAM(Points)
%   NeLCS = NEDLOADLCS_BEAM(Points,phi_y,phi_z)
%   determines the values of the interpolation functions for a distributed 
%   load in the specified points (LCS). These are used to compute the 
%   displacements that occur due to the distributed loads if all nodes are 
%   fixed.
%
%   NeLCS = NEDLOADLCS_BEAM(Points,[],[],a,b,L)
%   NeLCS = NEDLOADLCS_BEAM(Points,phi_y,phi_z,a,b,L)
%   determines the values of the interpolation functions for a partial
%   distributed load in the specified points (LCS). The load starts at a 
%   distance 'a' and ends at distance 'b' from the first node of the element:
%       1) from 0 to a: no DLoad,
%       2) from a to b: DLoad,
%       3) from b to L: no DLoad.
%
%   [NeLCS,dNeLCSdx] = NEDLOADLCS_BEAM(Points,[],[],a,b,L,dadx,dbdx,dLdx)
%   [NeLCS,dNeLCSdx] = NEDLOADLCS_BEAM(Points,phi_y,phi_z,a,b,L,dadx,dbdx,dLdx)
%   additionally computes the derivatives of the interpolation functions of
%   a partial distributed load with respect to the design variables x.
%   Note: the derivatives of the interpolation functions for a distributed
%   load on the complete element are zero.
%
%   Points     Points in the local coordinate system (1 * nPoints) 
%   phi_y      Shear deformation constant in y dir   (scalar)
%   phi_z      Shear deformation constant in z dir   (scalar)
%   a          Local starting point for the DLoad    (scalar)
%   b          Local ending point for the DLoad      (scalar)
%   L          Element length                        (scalar)
%   dadx       Local starting point derivatives      (1 * nVar)
%   dbdx       Local ending point derivatives        (1 * nVar)
%   dLdx       Element length derivatives            (1 * nVar)
%   NeDLoad    Interpolated values                   (nPoints * 6)
%   dNeDLoaddx Interpolated values derivatives       (nPoints * 6 *nVar)
%
%   See also DISP_BEAM, NELCS_BEAM.

% David Dooms, Wouter Dillen
% September 2008, December 2016

if nargin<2, phi_y = []; end
if nargin<3, phi_z = []; end
if nargin<4, a = []; end
if nargin<5, b = []; end
if nargin<6, L = []; end
if nargin<7, dadx = []; end
if nargin<8, dbdx = []; end
if nargin<9, dLdx = []; end

nVar = 0;
if nargout>1 && (~isempty(dadx) || ~isempty(dbdx) || ~isempty(dLdx))
    nVar = max([length(dadx),length(dbdx),length(dLdx)]);
end

if isempty(phi_y), phi_y = 0; end
if isempty(phi_z), phi_z = 0; end
if isempty(a), a = 0; end
if isempty(b), b = L; end

Points = Points(:).';
PartialDLoad = (~isempty(L) && (a~=0 || b~=L));


%% DLoad on the complete element %%
if not(PartialDLoad)

    if phi_y==0 && phi_z==0
        % Interpolation functions for a Bernoulli beam
        A=[ 0   0  0   0  0  0;
            -1  5  -7  3  0  0;
            -1  5  -7  3  0  0;
            0   0  0   0  0  0;
            1   0  -3  2  0  0;
            1   0  -3  2  0  0];
    else
        % Interpolation functions for a Timoshenko beam
        A=[ 0  0  0  0  0  0;
            -1  5 (5*phi_y^2-15*phi_y-21)/3/(1+phi_y)  (-10*phi_y^2-5*phi_y+6)/2/(1+phi_y)  phi_y*(20*phi_y+21)/6/(1+phi_y)  0;
            -1  5 (5*phi_z^2-15*phi_z-21)/3/(1+phi_z)  (-10*phi_z^2-5*phi_z+6)/2/(1+phi_z)  phi_z*(20*phi_z+21)/6/(1+phi_z)  0;
            0  0  0  0  0  0;
            1  0 (-5*phi_y^2-15*phi_y-9)/3/(1+phi_y)  (5*phi_y+4)/2/(1+phi_y)              phi_y*(10*phi_y+9)/6/(1+phi_y)   0;
            1  0 (-5*phi_z^2-15*phi_z-9)/3/(1+phi_z)  (5*phi_z+4)/2/(1+phi_z)              phi_z*(10*phi_z+9)/6/(1+phi_z)   0];
    end
    
    % NeDLoad
    NeDLoad=zeros(length(Points),6);
    for k=1:6
        NeDLoad(:,k)=polyval(A(k,:),Points);
    end
    
    % dNeDLoaddx
    if nargout>1
        dNeDLoaddx=zeros([size(NeDLoad),nVar]);
    end
    
end

%% DLoad on part of the element %%
if PartialDLoad
    
    if phi_y==0 && phi_z==0
        % Interpolation functions for a Bernoulli beam
        
        A1=[ 0  0  0  0  0  0;
             0  0  2*(4*a^4-b^4)-2*a^2*b^2-5*L*a*(3*a^2-2*L^2)+5*L*b*(b^2-2*L^2)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)  -(-10*L*a^2*(3*a-2*L)+10*L*b^2*(b-L)+3*(4*a^4-b^4)-3*a*b*(a^2+b^2)-10*L^2*a*b-3*a^2*b^2+10*L*a*b*(a+b))  0  0;
             0  0  2*(4*a^4-b^4)-2*a^2*b^2-5*L*a*(3*a^2-2*L^2)+5*L*b*(b^2-2*L^2)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)  -(-10*L*a^2*(3*a-2*L)+10*L*b^2*(b-L)+3*(4*a^4-b^4)-3*a*b*(a^2+b^2)-10*L^2*a*b-3*a^2*b^2+10*L*a*b*(a+b))  0  0;
             0  0  0  0  0  0;
             0  0  2*(a^4-4*b^4)+2*a^2*b^2-5*L*a*(a^2-2*L^2)+5*L*b*(3*b^2-2*L^2)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)  -(-10*L*a^2*(a-L)+10*L*b^2*(3*b-2*L)+3*(a^4-4*b^4)+3*a*b*(a^2+b^2)+10*L^2*a*b+3*a^2*b^2-10*L*a*b*(a+b))  0  0;
             0  0  2*(a^4-4*b^4)+2*a^2*b^2-5*L*a*(a^2-2*L^2)+5*L*b*(3*b^2-2*L^2)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)  -(-10*L*a^2*(a-L)+10*L*b^2*(3*b-2*L)+3*(a^4-4*b^4)+3*a*b*(a^2+b^2)+10*L^2*a*b+3*a^2*b^2-10*L*a*b*(a+b))  0  0];
        A2=[ 0  0  0  0  0  0;
             -L^5/(b-a)  5*L^4*b/(b-a)  (-2*(4*a^5+b^5)-10*L^3*b^2+5*L*(3*a^4+b^4)+10*a^4*b-20*L*a^3*b)/(b-a)  -(10*L*(3*a^4+b^4)-10*L^2*b^3-3*(4*a^5+b^5)+15*a^4*b-40*L*a^3*b)/(b-a)   5*L*(3*a^4-4*a^3*b)/(b-a)  (-4*a^5+5*a^4*b)/(b-a);
             -L^5/(b-a)  5*L^4*b/(b-a)  (-2*(4*a^5+b^5)-10*L^3*b^2+5*L*(3*a^4+b^4)+10*a^4*b-20*L*a^3*b)/(b-a)  -(10*L*(3*a^4+b^4)-10*L^2*b^3-3*(4*a^5+b^5)+15*a^4*b-40*L*a^3*b)/(b-a)   5*L*(3*a^4-4*a^3*b)/(b-a)  (-4*a^5+5*a^4*b)/(b-a);
             0  0  0  0  0  0;
             L^5/(b-a)  -5*L^4*a/(b-a)  (-2*(a^5+4*b^5)-10*L^3*b^2+5*L*(a^4+3*b^4)+10*a*b^4-20*L*a*b^3+20*L^3*a*b)/(b-a)  -(10*L*(a^4+3*b^4)-20*L^2*b^3-3*(a^5+4*b^5)+15*a*b^4-40*L*a*b^3+30*L^2*a*b^2)/(b-a)  5*L*a^4/(b-a)  -a^5/(b-a);
             L^5/(b-a)  -5*L^4*a/(b-a)  (-2*(a^5+4*b^5)-10*L^3*b^2+5*L*(a^4+3*b^4)+10*a*b^4-20*L*a*b^3+20*L^3*a*b)/(b-a)  -(10*L*(a^4+3*b^4)-20*L^2*b^3-3*(a^5+4*b^5)+15*a*b^4-40*L*a*b^3+30*L^2*a*b^2)/(b-a)  5*L*a^4/(b-a)  -a^5/(b-a)];
        A3=[ 0  0  0  0  0  0;
             0  0  2*(4*a^4-b^4)-2*a^2*b^2-5*L*(3*a^3-b^3)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)  -(-10*L*(3*a^3-b^3)+3*(4*a^4-b^4)-3*a*b*(a^2+b^2)-3*a^2*b^2+10*L*a*b*(a+b))  5*L*(b^3+a*b^2+a^2*b-3*a^3)  (4*a^4-b^4)-a^2*b^2-a*b*(a^2+b^2);
             0  0  2*(4*a^4-b^4)-2*a^2*b^2-5*L*(3*a^3-b^3)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)  -(-10*L*(3*a^3-b^3)+3*(4*a^4-b^4)-3*a*b*(a^2+b^2)-3*a^2*b^2+10*L*a*b*(a+b))  5*L*(b^3+a*b^2+a^2*b-3*a^3)  (4*a^4-b^4)-a^2*b^2-a*b*(a^2+b^2);
             0  0  0  0  0  0;
             0  0  2*(a^4-4*b^4)+2*a^2*b^2-5*L*(a^3-3*b^3)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)  -(-10*L*(a^3-3*b^3)+3*(a^4-4*b^4)+3*a*b*(a^2+b^2)+3*a^2*b^2-10*L*a*b*(a+b))  5*L*(3*b^3-a*b^2-a^2*b-a^3)  (a^4-4*b^4)+a^2*b^2+a*b*(a^2+b^2);
             0  0  2*(a^4-4*b^4)+2*a^2*b^2-5*L*(a^3-3*b^3)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)  -(-10*L*(a^3-3*b^3)+3*(a^4-4*b^4)+3*a*b*(a^2+b^2)+3*a^2*b^2-10*L*a*b*(a+b))  5*L*(3*b^3-a*b^2-a^2*b-a^3)  (a^4-4*b^4)+a^2*b^2+a*b*(a^2+b^2)];
      
        dA1dx = zeros(6,6,nVar);
        dA2dx = zeros(6,6,nVar);
        dA3dx = zeros(6,6,nVar);
        for n=1:nVar
            dA1dx(:,:,n)=[ 0  0  0  0  0  0;
                           0  0  2*(16*a^3*dadx(n)-4*b^3*dbdx(n))-4*a*dadx(n)*b^2-4*a^2*b*dbdx(n)-5*dLdx(n)*a*(3*a^2-2*L^2)-5*L*dadx(n)*(3*a^2-2*L^2)-5*L*a*(6*a*dadx(n)-4*L*dLdx(n))+5*dLdx(n)*b*(b^2-2*L^2)+5*L*dbdx(n)*(b^2-2*L^2)+5*L*b*(2*b*dbdx(n)-4*L*dLdx(n))-2*dadx(n)*b*(a^2+b^2)-2*a*dbdx(n)*(a^2+b^2)-2*a*b*(2*a*dadx(n)+2*b*dbdx(n))+5*dLdx(n)*a*b*(a+b)+5*L*dadx(n)*b*(a+b)+5*L*a*dbdx(n)*(a+b)+5*L*a*b*(dadx(n)+dbdx(n))  -(-10*dLdx(n)*a^2*(3*a-2*L)-20*L*a*dadx(n)*(3*a-2*L)-10*L*a^2*(3*dadx(n)-2*dLdx(n))+10*dLdx(n)*b^2*(b-L)+20*L*b*dbdx(n)*(b-L)+10*L*b^2*(dbdx(n)-dLdx(n))+3*(16*a^3*dadx(n)-4*b^3*dbdx(n))-3*dadx(n)*b*(a^2+b^2)-3*a*dbdx(n)*(a^2+b^2)-3*a*b*(2*a*dadx(n)+2*b*dbdx(n))-20*L*dLdx(n)*a*b-10*L^2*dadx(n)*b-10*L^2*a*dbdx(n)-6*a*dadx(n)*b^2-6*a^2*b*dbdx(n)+10*dLdx(n)*a*b*(a+b)+10*L*dadx(n)*b*(a+b)+10*L*a*dbdx(n)*(a+b)+10*L*a*b*(dadx(n)+dbdx(n)))  0  0;
                           0  0  2*(16*a^3*dadx(n)-4*b^3*dbdx(n))-4*a*dadx(n)*b^2-4*a^2*b*dbdx(n)-5*dLdx(n)*a*(3*a^2-2*L^2)-5*L*dadx(n)*(3*a^2-2*L^2)-5*L*a*(6*a*dadx(n)-4*L*dLdx(n))+5*dLdx(n)*b*(b^2-2*L^2)+5*L*dbdx(n)*(b^2-2*L^2)+5*L*b*(2*b*dbdx(n)-4*L*dLdx(n))-2*dadx(n)*b*(a^2+b^2)-2*a*dbdx(n)*(a^2+b^2)-2*a*b*(2*a*dadx(n)+2*b*dbdx(n))+5*dLdx(n)*a*b*(a+b)+5*L*dadx(n)*b*(a+b)+5*L*a*dbdx(n)*(a+b)+5*L*a*b*(dadx(n)+dbdx(n))  -(-10*dLdx(n)*a^2*(3*a-2*L)-20*L*a*dadx(n)*(3*a-2*L)-10*L*a^2*(3*dadx(n)-2*dLdx(n))+10*dLdx(n)*b^2*(b-L)+20*L*b*dbdx(n)*(b-L)+10*L*b^2*(dbdx(n)-dLdx(n))+3*(16*a^3*dadx(n)-4*b^3*dbdx(n))-3*dadx(n)*b*(a^2+b^2)-3*a*dbdx(n)*(a^2+b^2)-3*a*b*(2*a*dadx(n)+2*b*dbdx(n))-20*L*dLdx(n)*a*b-10*L^2*dadx(n)*b-10*L^2*a*dbdx(n)-6*a*dadx(n)*b^2-6*a^2*b*dbdx(n)+10*dLdx(n)*a*b*(a+b)+10*L*dadx(n)*b*(a+b)+10*L*a*dbdx(n)*(a+b)+10*L*a*b*(dadx(n)+dbdx(n)))  0  0;
                           0  0  0  0  0  0;
                           0  0  2*(4*a^3*dadx(n)-16*b^3*dbdx(n))+4*a*dadx(n)*b^2+4*a^2*b*dbdx(n)-5*dLdx(n)*a*(a^2-2*L^2)-5*L*dadx(n)*(a^2-2*L^2)-5*L*a*(2*a*dadx(n)-4*L*dLdx(n))+5*dLdx(n)*b*(3*b^2-2*L^2)+5*L*dbdx(n)*(3*b^2-2*L^2)+5*L*b*(6*b*dbdx(n)-4*L*dLdx(n))+2*dadx(n)*b*(a^2+b^2)+2*a*dbdx(n)*(a^2+b^2)+2*a*b*(2*a*dadx(n)+2*b*dbdx(n))-5*dLdx(n)*a*b*(a+b)-5*L*dadx(n)*b*(a+b)-5*L*a*dbdx(n)*(a+b)-5*L*a*b*(dadx(n)+dbdx(n))  -(-10*dLdx(n)*a^2*(a-L)-20*L*a*dadx(n)*(a-L)-10*L*a^2*(dadx(n)-dLdx(n))+10*dLdx(n)*b^2*(3*b-2*L)+20*L*b*dbdx(n)*(3*b-2*L)+10*L*b^2*(3*dbdx(n)-2*dLdx(n))+3*(4*a^3*dadx(n)-16*b^3*dbdx(n))+3*dadx(n)*b*(a^2+b^2)+3*a*dbdx(n)*(a^2+b^2)+3*a*b*(2*a*dadx(n)+2*b*dbdx(n))+20*L*dLdx(n)*a*b+10*L^2*dadx(n)*b+10*L^2*a*dbdx(n)+6*a*dadx(n)*b^2+6*a^2*b*dbdx(n)-10*dLdx(n)*a*b*(a+b)-10*L*dadx(n)*b*(a+b)-10*L*a*dbdx(n)*(a+b)-10*L*a*b*(dadx(n)+dbdx(n)))  0  0;
                           0  0  2*(4*a^3*dadx(n)-16*b^3*dbdx(n))+4*a*dadx(n)*b^2+4*a^2*b*dbdx(n)-5*dLdx(n)*a*(a^2-2*L^2)-5*L*dadx(n)*(a^2-2*L^2)-5*L*a*(2*a*dadx(n)-4*L*dLdx(n))+5*dLdx(n)*b*(3*b^2-2*L^2)+5*L*dbdx(n)*(3*b^2-2*L^2)+5*L*b*(6*b*dbdx(n)-4*L*dLdx(n))+2*dadx(n)*b*(a^2+b^2)+2*a*dbdx(n)*(a^2+b^2)+2*a*b*(2*a*dadx(n)+2*b*dbdx(n))-5*dLdx(n)*a*b*(a+b)-5*L*dadx(n)*b*(a+b)-5*L*a*dbdx(n)*(a+b)-5*L*a*b*(dadx(n)+dbdx(n))  -(-10*dLdx(n)*a^2*(a-L)-20*L*a*dadx(n)*(a-L)-10*L*a^2*(dadx(n)-dLdx(n))+10*dLdx(n)*b^2*(3*b-2*L)+20*L*b*dbdx(n)*(3*b-2*L)+10*L*b^2*(3*dbdx(n)-2*dLdx(n))+3*(4*a^3*dadx(n)-16*b^3*dbdx(n))+3*dadx(n)*b*(a^2+b^2)+3*a*dbdx(n)*(a^2+b^2)+3*a*b*(2*a*dadx(n)+2*b*dbdx(n))+20*L*dLdx(n)*a*b+10*L^2*dadx(n)*b+10*L^2*a*dbdx(n)+6*a*dadx(n)*b^2+6*a^2*b*dbdx(n)-10*dLdx(n)*a*b*(a+b)-10*L*dadx(n)*b*(a+b)-10*L*a*dbdx(n)*(a+b)-10*L*a*b*(dadx(n)+dbdx(n)))  0  0];
            dA2dx(:,:,n)=[ 0  0  0  0  0  0;
                           -5*L^4*dLdx(n)/(b-a)-L^5*(dadx(n)-dbdx(n))/(b-a)^2  20*L^3*dLdx(n)*b/(b-a)+5*L^4*dbdx(n)/(b-a)+5*L^4*b*(dadx(n)-dbdx(n))/(b-a)^2  (-2*(20*a^4*dadx(n)+5*b^4*dbdx(n))-30*L^2*dLdx(n)*b^2-20*L^3*b*dbdx(n)+5*dLdx(n)*(3*a^4+b^4)+5*L*(12*a^3*dadx(n)+4*b^3*dbdx(n))+40*a^3*dadx(n)*b+10*a^4*dbdx(n)-20*dLdx(n)*a^3*b-60*L*a^2*dadx(n)*b-20*L*a^3*dbdx(n))/(b-a)+(-2*(4*a^5+b^5)-10*L^3*b^2+5*L*(3*a^4+b^4)+10*a^4*b-20*L*a^3*b)*(dadx(n)-dbdx(n))/(b-a)^2  -(10*dLdx(n)*(3*a^4+b^4)+10*L*(12*a^3*dadx(n)+4*b^3*dbdx(n))-20*L*dLdx(n)*b^3-30*L^2*b^2*dbdx(n)-3*(20*a^4*dadx(n)+5*b^4*dbdx(n))+60*a^3*dadx(n)*b+15*a^4*dbdx(n)-40*dLdx(n)*a^3*b-120*L*a^2*dadx(n)*b-40*L*a^3*dbdx(n))/(b-a)-(10*L*(3*a^4+b^4)-10*L^2*b^3-3*(4*a^5+b^5)+15*a^4*b-40*L*a^3*b)*(dadx(n)-dbdx(n))/(b-a)^2  (5*dLdx(n)*(3*a^4-4*a^3*b)+5*L*(12*a^3*dadx(n)-12*a^2*dadx(n)*b-4*a^3*dbdx(n)))/(b-a)+5*L*(3*a^4-4*a^3*b)*(dadx(n)-dbdx(n))/(b-a)^2  (-20*a^4*dadx(n)+20*a^3*dadx(n)*b+5*a^4*dbdx(n))/(b-a)+(-4*a^5+5*a^4*b)*(dadx(n)-dbdx(n))/(b-a)^2;
                           -5*L^4*dLdx(n)/(b-a)-L^5*(dadx(n)-dbdx(n))/(b-a)^2  20*L^3*dLdx(n)*b/(b-a)+5*L^4*dbdx(n)/(b-a)+5*L^4*b*(dadx(n)-dbdx(n))/(b-a)^2  (-2*(20*a^4*dadx(n)+5*b^4*dbdx(n))-30*L^2*dLdx(n)*b^2-20*L^3*b*dbdx(n)+5*dLdx(n)*(3*a^4+b^4)+5*L*(12*a^3*dadx(n)+4*b^3*dbdx(n))+40*a^3*dadx(n)*b+10*a^4*dbdx(n)-20*dLdx(n)*a^3*b-60*L*a^2*dadx(n)*b-20*L*a^3*dbdx(n))/(b-a)+(-2*(4*a^5+b^5)-10*L^3*b^2+5*L*(3*a^4+b^4)+10*a^4*b-20*L*a^3*b)*(dadx(n)-dbdx(n))/(b-a)^2  -(10*dLdx(n)*(3*a^4+b^4)+10*L*(12*a^3*dadx(n)+4*b^3*dbdx(n))-20*L*dLdx(n)*b^3-30*L^2*b^2*dbdx(n)-3*(20*a^4*dadx(n)+5*b^4*dbdx(n))+60*a^3*dadx(n)*b+15*a^4*dbdx(n)-40*dLdx(n)*a^3*b-120*L*a^2*dadx(n)*b-40*L*a^3*dbdx(n))/(b-a)-(10*L*(3*a^4+b^4)-10*L^2*b^3-3*(4*a^5+b^5)+15*a^4*b-40*L*a^3*b)*(dadx(n)-dbdx(n))/(b-a)^2  (5*dLdx(n)*(3*a^4-4*a^3*b)+5*L*(12*a^3*dadx(n)-12*a^2*dadx(n)*b-4*a^3*dbdx(n)))/(b-a)+5*L*(3*a^4-4*a^3*b)*(dadx(n)-dbdx(n))/(b-a)^2  (-20*a^4*dadx(n)+20*a^3*dadx(n)*b+5*a^4*dbdx(n))/(b-a)+(-4*a^5+5*a^4*b)*(dadx(n)-dbdx(n))/(b-a)^2;
                           0  0  0  0  0  0;
                           5*L^4*dLdx(n)/(b-a)+L^5*(dadx(n)-dbdx(n))/(b-a)^2  -20*L^3*dLdx(n)*a/(b-a)-5*L^4*dadx(n)/(b-a)-5*L^4*a*(dadx(n)-dbdx(n))/(b-a)^2  (-2*(5*a^4*dadx(n)+20*b^4*dbdx(n))-30*L^2*dLdx(n)*b^2-20*L^3*b*dbdx(n)+5*dLdx(n)*(a^4+3*b^4)+5*L*(4*a^3*dadx(n)+12*b^3*dbdx(n))+10*dadx(n)*b^4+40*a*b^3*dbdx(n)-20*dLdx(n)*a*b^3-20*L*dadx(n)*b^3-60*L*a*b^2*dbdx(n)+60*L^2*dLdx(n)*a*b+20*L^3*dadx(n)*b+20*L^3*a*dbdx(n))/(b-a)+(-2*(a^5+4*b^5)-10*L^3*b^2+5*L*(a^4+3*b^4)+10*a*b^4-20*L*a*b^3+20*L^3*a*b)*(dadx(n)-dbdx(n))/(b-a)^2  -(10*dLdx(n)*(a^4+3*b^4)+10*L*(4*a^3*dadx(n)+12*b^3*dbdx(n))-40*L*dLdx(n)*b^3-60*L^2*b^2*dbdx(n)-3*(5*a^4*dadx(n)+20*b^4*dbdx(n))+15*dadx(n)*b^4+60*a*b^3*dbdx(n)-40*dLdx(n)*a*b^3-40*L*dadx(n)*b^3-120*L*a*b^2*dbdx(n)+60*L*dLdx(n)*a*b^2+30*L^2*dadx(n)*b^2+60*L^2*a*b*dbdx(n))/(b-a)-(10*L*(a^4+3*b^4)-20*L^2*b^3-3*(a^5+4*b^5)+15*a*b^4-40*L*a*b^3+30*L^2*a*b^2)*(dadx(n)-dbdx(n))/(b-a)^2  (5*dLdx(n)*a^4+20*L*a^3*dadx(n))/(b-a)+5*L*a^4*(dadx(n)-dbdx(n))/(b-a)^2  (-5*a^4*dadx(n))/(b-a)+(-a^5)*(dadx(n)-dbdx(n))/(b-a)^2;
                           5*L^4*dLdx(n)/(b-a)+L^5*(dadx(n)-dbdx(n))/(b-a)^2  -20*L^3*dLdx(n)*a/(b-a)-5*L^4*dadx(n)/(b-a)-5*L^4*a*(dadx(n)-dbdx(n))/(b-a)^2  (-2*(5*a^4*dadx(n)+20*b^4*dbdx(n))-30*L^2*dLdx(n)*b^2-20*L^3*b*dbdx(n)+5*dLdx(n)*(a^4+3*b^4)+5*L*(4*a^3*dadx(n)+12*b^3*dbdx(n))+10*dadx(n)*b^4+40*a*b^3*dbdx(n)-20*dLdx(n)*a*b^3-20*L*dadx(n)*b^3-60*L*a*b^2*dbdx(n)+60*L^2*dLdx(n)*a*b+20*L^3*dadx(n)*b+20*L^3*a*dbdx(n))/(b-a)+(-2*(a^5+4*b^5)-10*L^3*b^2+5*L*(a^4+3*b^4)+10*a*b^4-20*L*a*b^3+20*L^3*a*b)*(dadx(n)-dbdx(n))/(b-a)^2  -(10*dLdx(n)*(a^4+3*b^4)+10*L*(4*a^3*dadx(n)+12*b^3*dbdx(n))-40*L*dLdx(n)*b^3-60*L^2*b^2*dbdx(n)-3*(5*a^4*dadx(n)+20*b^4*dbdx(n))+15*dadx(n)*b^4+60*a*b^3*dbdx(n)-40*dLdx(n)*a*b^3-40*L*dadx(n)*b^3-120*L*a*b^2*dbdx(n)+60*L*dLdx(n)*a*b^2+30*L^2*dadx(n)*b^2+60*L^2*a*b*dbdx(n))/(b-a)-(10*L*(a^4+3*b^4)-20*L^2*b^3-3*(a^5+4*b^5)+15*a*b^4-40*L*a*b^3+30*L^2*a*b^2)*(dadx(n)-dbdx(n))/(b-a)^2  (5*dLdx(n)*a^4+20*L*a^3*dadx(n))/(b-a)+5*L*a^4*(dadx(n)-dbdx(n))/(b-a)^2  (-5*a^4*dadx(n))/(b-a)+(-a^5)*(dadx(n)-dbdx(n))/(b-a)^2];
            dA3dx(:,:,n)=[ 0  0  0  0  0  0;
                           0  0  2*(16*a^3*dadx(n)-4*b^3*dbdx(n))-4*a*dadx(n)*b^2-4*a^2*b*dbdx(n)-5*dLdx(n)*(3*a^3-b^3)-5*L*(9*a^2*dadx(n)-3*b^2*dbdx(n))-2*dadx(n)*b*(a^2+b^2)-2*a*dbdx(n)*(a^2+b^2)-2*a*b*(2*a*dadx(n)+2*b*dbdx(n))+5*dLdx(n)*a*b*(a+b)+5*L*dadx(n)*b*(a+b)+5*L*a*dbdx(n)*(a+b)+5*L*a*b*(dadx(n)+dbdx(n))  -(-10*dLdx(n)*(3*a^3-b^3)-10*L*(9*a^2*dadx(n)-3*b^2*dbdx(n))+3*(16*a^3*dadx(n)-4*b^3*dbdx(n))-3*dadx(n)*b*(a^2+b^2)-3*a*dbdx(n)*(a^2+b^2)-3*a*b*(2*a*dadx(n)+2*b*dbdx(n))-6*a*dadx(n)*b^2-6*a^2*b*dbdx(n)+10*dLdx(n)*a*b*(a+b)+10*L*dadx(n)*b*(a+b)+10*L*a*dbdx(n)*(a+b)+10*L*a*b*(dadx(n)+dbdx(n)))  5*dLdx(n)*(b^3+a*b^2+a^2*b-3*a^3)+5*L*(3*b^2*dbdx(n)+dadx(n)*b^2+2*a*b*(dadx(n)+dbdx(n))+a^2*dbdx(n)-9*a^2*dadx(n))  (16*a^3*dadx(n)-4*b^3*dbdx(n))-2*a*dadx(n)*b^2-2*a^2*b*dbdx(n)-dadx(n)*b*(a^2+b^2)-a*dbdx(n)*(a^2+b^2)-a*b*(2*a*dadx(n)+2*b*dbdx(n));
                           0  0  2*(16*a^3*dadx(n)-4*b^3*dbdx(n))-4*a*dadx(n)*b^2-4*a^2*b*dbdx(n)-5*dLdx(n)*(3*a^3-b^3)-5*L*(9*a^2*dadx(n)-3*b^2*dbdx(n))-2*dadx(n)*b*(a^2+b^2)-2*a*dbdx(n)*(a^2+b^2)-2*a*b*(2*a*dadx(n)+2*b*dbdx(n))+5*dLdx(n)*a*b*(a+b)+5*L*dadx(n)*b*(a+b)+5*L*a*dbdx(n)*(a+b)+5*L*a*b*(dadx(n)+dbdx(n))  -(-10*dLdx(n)*(3*a^3-b^3)-10*L*(9*a^2*dadx(n)-3*b^2*dbdx(n))+3*(16*a^3*dadx(n)-4*b^3*dbdx(n))-3*dadx(n)*b*(a^2+b^2)-3*a*dbdx(n)*(a^2+b^2)-3*a*b*(2*a*dadx(n)+2*b*dbdx(n))-6*a*dadx(n)*b^2-6*a^2*b*dbdx(n)+10*dLdx(n)*a*b*(a+b)+10*L*dadx(n)*b*(a+b)+10*L*a*dbdx(n)*(a+b)+10*L*a*b*(dadx(n)+dbdx(n)))  5*dLdx(n)*(b^3+a*b^2+a^2*b-3*a^3)+5*L*(3*b^2*dbdx(n)+dadx(n)*b^2+2*a*b*(dadx(n)+dbdx(n))+a^2*dbdx(n)-9*a^2*dadx(n))  (16*a^3*dadx(n)-4*b^3*dbdx(n))-2*a*dadx(n)*b^2-2*a^2*b*dbdx(n)-dadx(n)*b*(a^2+b^2)-a*dbdx(n)*(a^2+b^2)-a*b*(2*a*dadx(n)+2*b*dbdx(n));
                           0  0  0  0  0  0;
                           0  0  2*(4*a^3*dadx(n)-16*b^3*dbdx(n))+4*a*dadx(n)*b^2+4*a^2*b*dbdx(n)-5*dLdx(n)*(a^3-3*b^3)-5*L*(3*a^2*dadx(n)-9*b^2*dbdx(n))+2*dadx(n)*b*(a^2+b^2)+2*a*dbdx(n)*(a^2+b^2)+2*a*b*(2*a*dadx(n)+2*b*dbdx(n))-5*dLdx(n)*a*b*(a+b)-5*L*dadx(n)*b*(a+b)-5*L*a*dbdx(n)*(a+b)-5*L*a*b*(dadx(n)+dbdx(n))  -(-10*dLdx(n)*(a^3-3*b^3)-10*L*(3*a^2*dadx(n)-9*b^2*dbdx(n))+3*(4*a^3*dadx(n)-16*b^3*dbdx(n))+3*dadx(n)*b*(a^2+b^2)+3*a*dbdx(n)*(a^2+b^2)+3*a*b*(2*a*dadx(n)+2*b*dbdx(n))+6*a*dadx(n)*b^2+6*a^2*b*dbdx(n)-10*dLdx(n)*a*b*(a+b)-10*L*dadx(n)*b*(a+b)-10*L*a*dbdx(n)*(a+b)-10*L*a*b*(dadx(n)+dbdx(n)))  5*dLdx(n)*(3*b^3-a*b^2-a^2*b-a^3)+5*L*(9*b^2*dbdx(n)-dadx(n)*b^2-2*a*b*(dadx(n)+dbdx(n))-a^2*dbdx(n)-3*a^2*dadx(n))  (4*a^3*dadx(n)-16*b^3*dbdx(n))+2*a*dadx(n)*b^2+2*a^2*b*dbdx(n)+dadx(n)*b*(a^2+b^2)+a*dbdx(n)*(a^2+b^2)+a*b*(2*a*dadx(n)+2*b*dbdx(n));
                           0  0  2*(4*a^3*dadx(n)-16*b^3*dbdx(n))+4*a*dadx(n)*b^2+4*a^2*b*dbdx(n)-5*dLdx(n)*(a^3-3*b^3)-5*L*(3*a^2*dadx(n)-9*b^2*dbdx(n))+2*dadx(n)*b*(a^2+b^2)+2*a*dbdx(n)*(a^2+b^2)+2*a*b*(2*a*dadx(n)+2*b*dbdx(n))-5*dLdx(n)*a*b*(a+b)-5*L*dadx(n)*b*(a+b)-5*L*a*dbdx(n)*(a+b)-5*L*a*b*(dadx(n)+dbdx(n))  -(-10*dLdx(n)*(a^3-3*b^3)-10*L*(3*a^2*dadx(n)-9*b^2*dbdx(n))+3*(4*a^3*dadx(n)-16*b^3*dbdx(n))+3*dadx(n)*b*(a^2+b^2)+3*a*dbdx(n)*(a^2+b^2)+3*a*b*(2*a*dadx(n)+2*b*dbdx(n))+6*a*dadx(n)*b^2+6*a^2*b*dbdx(n)-10*dLdx(n)*a*b*(a+b)-10*L*dadx(n)*b*(a+b)-10*L*a*dbdx(n)*(a+b)-10*L*a*b*(dadx(n)+dbdx(n)))  5*dLdx(n)*(3*b^3-a*b^2-a^2*b-a^3)+5*L*(9*b^2*dbdx(n)-dadx(n)*b^2-2*a*b*(dadx(n)+dbdx(n))-a^2*dbdx(n)-3*a^2*dadx(n))  (4*a^3*dadx(n)-16*b^3*dbdx(n))+2*a*dadx(n)*b^2+2*a^2*b*dbdx(n)+dadx(n)*b*(a^2+b^2)+a*dbdx(n)*(a^2+b^2)+a*b*(2*a*dadx(n)+2*b*dbdx(n))];
        end
    
    else
        % Interpolation functions for a Timoshenko beam
   
        A1=[ 0  0  0  0  0  0;
             0  0  (phi_y*L^2/3*(-10*(2*a^2-b^2)+30*L*(a-b)+10*a*b)+(2*(4*a^4-b^4)-2*a^2*b^2-5*L*a*(3*a^2-2*L^2)+5*L*b*(b^2-2*L^2)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)))/(1+phi_y)  -(phi_y*L/4*(-10*(3*a^3-b^3)+20*L*(2*a^2-b^2)+10*a*b*(a+b)-20*L*a*b)+(-10*L*a^2*(3*a-2*L)+10*L*b^2*(b-L)+3*(4*a^4-b^4)-3*a*b*(a^2+b^2)-10*L^2*a*b-3*a^2*b^2+10*L*a*b*(a+b)))/(1+phi_y)  -((phi_y*L)^2/6*(-10*(2*a^2-b^2)+30*L*(a-b)+10*a*b)+phi_y/2*(2*(4*a^4-b^4)-2*a^2*b^2-5*L*a*(3*a^2-2*L^2)+5*L*b*(b^2-2*L^2)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)))/(1+phi_y)  0;
             0  0  (phi_z*L^2/3*(-10*(2*a^2-b^2)+30*L*(a-b)+10*a*b)+(2*(4*a^4-b^4)-2*a^2*b^2-5*L*a*(3*a^2-2*L^2)+5*L*b*(b^2-2*L^2)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)))/(1+phi_z)  -(phi_z*L/4*(-10*(3*a^3-b^3)+20*L*(2*a^2-b^2)+10*a*b*(a+b)-20*L*a*b)+(-10*L*a^2*(3*a-2*L)+10*L*b^2*(b-L)+3*(4*a^4-b^4)-3*a*b*(a^2+b^2)-10*L^2*a*b-3*a^2*b^2+10*L*a*b*(a+b)))/(1+phi_z)  -((phi_z*L)^2/6*(-10*(2*a^2-b^2)+30*L*(a-b)+10*a*b)+phi_z/2*(2*(4*a^4-b^4)-2*a^2*b^2-5*L*a*(3*a^2-2*L^2)+5*L*b*(b^2-2*L^2)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)))/(1+phi_z)  0;
             0  0  0  0  0  0;
             0  0  (phi_y*L^2/3*(-10*(a^2-2*b^2)+30*L*(a-b)-10*a*b)+(2*(a^4-4*b^4)+2*a^2*b^2-5*L*a*(a^2-2*L^2)+5*L*b*(3*b^2-2*L^2)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)))/(1+phi_y)  -(phi_y*L/4*(-10*(a^3-3*b^3)+20*L*(a^2-2*b^2)-10*a*b*(a+b)+20*L*a*b)+(-10*L*a^2*(a-L)+10*L*b^2*(3*b-2*L)+3*(a^4-4*b^4)+3*a*b*(a^2+b^2)+10*L^2*a*b+3*a^2*b^2-10*L*a*b*(a+b)))/(1+phi_y)  -((phi_y*L)^2/6*(-10*(a^2-2*b^2)+30*L*(a-b)-10*a*b)+phi_y/2*(2*(a^4-4*b^4)+2*a^2*b^2-5*L*a*(a^2-2*L^2)+5*L*b*(3*b^2-2*L^2)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)))/(1+phi_y)  0;
             0  0  (phi_z*L^2/3*(-10*(a^2-2*b^2)+30*L*(a-b)-10*a*b)+(2*(a^4-4*b^4)+2*a^2*b^2-5*L*a*(a^2-2*L^2)+5*L*b*(3*b^2-2*L^2)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)))/(1+phi_z)  -(phi_z*L/4*(-10*(a^3-3*b^3)+20*L*(a^2-2*b^2)-10*a*b*(a+b)+20*L*a*b)+(-10*L*a^2*(a-L)+10*L*b^2*(3*b-2*L)+3*(a^4-4*b^4)+3*a*b*(a^2+b^2)+10*L^2*a*b+3*a^2*b^2-10*L*a*b*(a+b)))/(1+phi_z)  -((phi_z*L)^2/6*(-10*(a^2-2*b^2)+30*L*(a-b)-10*a*b)+phi_z/2*(2*(a^4-4*b^4)+2*a^2*b^2-5*L*a*(a^2-2*L^2)+5*L*b*(3*b^2-2*L^2)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)))/(1+phi_z)  0];       
        A2=[ 0  0  0  0  0  0;
             -L^5/(b-a)  5*L^4*b/(b-a)  (phi_y*L^2/3*(10*(2*a^3+b^3)-30*b*(L*b+a^2))+(-2*(4*a^5+b^5)-10*L^3*b^2+5*L*(3*a^4+b^4)+10*a^4*b-20*L*a^3*b))/(1+phi_y)/(b-a)+10/6*phi_y*L^5/(b-a)  -(phi_y*L/4*(10*(3*a^4+b^4)+20*L*(2*a^3-b^3)-40*a^3*b-60*L*a^2*b)+(10*L*(3*a^4+b^4)-10*L^2*b^3-3*(4*a^5+b^5)+15*a^4*b-40*L*a^3*b))/(1+phi_y)/(b-a)-5*phi_y*L^4*b/(b-a)  5*L*(3*a^4-4*a^3*b)/(b-a)-((phi_y*L)^2/6*(10*(2*a^3+b^3)-30*b*(L*b+a^2))+phi_y/2*(-2*(4*a^5+b^5)-10*L^3*b^2+5*L*(3*a^4+b^4)+10*a^4*b-20*L*a^3*b))/(1+phi_y)/(b-a)  (phi_y*L^2/6*(20*a^3-30*a^2*b)+(-4*a^5+5*a^4*b))/(b-a);
             -L^5/(b-a)  5*L^4*b/(b-a)  (phi_z*L^2/3*(10*(2*a^3+b^3)-30*b*(L*b+a^2))+(-2*(4*a^5+b^5)-10*L^3*b^2+5*L*(3*a^4+b^4)+10*a^4*b-20*L*a^3*b))/(1+phi_z)/(b-a)+10/6*phi_z*L^5/(b-a)  -(phi_z*L/4*(10*(3*a^4+b^4)+20*L*(2*a^3-b^3)-40*a^3*b-60*L*a^2*b)+(10*L*(3*a^4+b^4)-10*L^2*b^3-3*(4*a^5+b^5)+15*a^4*b-40*L*a^3*b))/(1+phi_z)/(b-a)-5*phi_z*L^4*b/(b-a)  5*L*(3*a^4-4*a^3*b)/(b-a)-((phi_z*L)^2/6*(10*(2*a^3+b^3)-30*b*(L*b+a^2))+phi_z/2*(-2*(4*a^5+b^5)-10*L^3*b^2+5*L*(3*a^4+b^4)+10*a^4*b-20*L*a^3*b))/(1+phi_z)/(b-a)  (phi_z*L^2/6*(20*a^3-30*a^2*b)+(-4*a^5+5*a^4*b))/(b-a);
             0  0  0  0  0  0;
             L^5/(b-a)  -5*L^4*a/(b-a)  (phi_y*L^2/3*(10*(a^3+2*b^3)-30*b^2*(a+L)+60*L*a*b)+(-2*(a^5+4*b^5)-10*L^3*b^2+5*L*(a^4+3*b^4)+10*a*b^4-20*L*a*b^3+20*L^3*a*b))/(1+phi_y)/(b-a)-10/6*phi_y*L^5/(b-a)  -(phi_y*L/4*(10*(a^4+3*b^4)+20*L*(a^3-2*b^3)-40*a*b^3+60*L*a*b^2)+(10*L*(a^4+3*b^4)-20*L^2*b^3-3*(a^5+4*b^5)+15*a*b^4-40*L*a*b^3+30*L^2*a*b^2))/(1+phi_y)/(b-a)+5*phi_y*L^4*a/(b-a)  5*L*a^4/(b-a)-((phi_y*L)^2/6*(10*(a^3+2*b^3)-30*b^2*(a+L)+60*L*a*b)+phi_y/2*(-2*(a^5+4*b^5)-10*L^3*b^2+5*L*(a^4+3*b^4)+10*a*b^4-20*L*a*b^3+20*L^3*a*b))/(1+phi_y)/(b-a)  (phi_y*L^2/6*10*a^3-a^5)/(b-a);
             L^5/(b-a)  -5*L^4*a/(b-a)  (phi_z*L^2/3*(10*(a^3+2*b^3)-30*b^2*(a+L)+60*L*a*b)+(-2*(a^5+4*b^5)-10*L^3*b^2+5*L*(a^4+3*b^4)+10*a*b^4-20*L*a*b^3+20*L^3*a*b))/(1+phi_z)/(b-a)-10/6*phi_z*L^5/(b-a)  -(phi_z*L/4*(10*(a^4+3*b^4)+20*L*(a^3-2*b^3)-40*a*b^3+60*L*a*b^2)+(10*L*(a^4+3*b^4)-20*L^2*b^3-3*(a^5+4*b^5)+15*a*b^4-40*L*a*b^3+30*L^2*a*b^2))/(1+phi_z)/(b-a)+5*phi_z*L^4*a/(b-a)  5*L*a^4/(b-a)-((phi_z*L)^2/6*(10*(a^3+2*b^3)-30*b^2*(a+L)+60*L*a*b)+phi_z/2*(-2*(a^5+4*b^5)-10*L^3*b^2+5*L*(a^4+3*b^4)+10*a*b^4-20*L*a*b^3+20*L^3*a*b))/(1+phi_z)/(b-a)  (phi_z*L^2/6*10*a^3-a^5)/(b-a)];
        A3=[ 0  0  0  0  0  0;    
             0  0  (phi_y*L^2/3*(-10*(2*a^2-b^2)+10*a*b)+(2*(4*a^4-b^4)-2*a^2*b^2-5*L*(3*a^3-b^3)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)))/(1+phi_y)  -(phi_y*L/4*(-10*(3*a^3-b^3)-20*L*(2*a^2-b^2)+10*a*b*(a+b)+20*L*a*b)+(-10*L*(3*a^3-b^3)+3*(4*a^4-b^4)-3*a*b*(a^2+b^2)-3*a^2*b^2+10*L*a*b*(a+b)))/(1+phi_y)  5*L*(b^3+a*b^2+a^2*b-3*a^3)-((phi_y*L)^2/6*(-10*(2*a^2-b^2)+10*a*b)+phi_y/2*(2*(4*a^4-b^4)-2*a^2*b^2-5*L*(3*a^3-b^3)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)))/(1+phi_y)  phi_y*L^2/6*(-10*(2*a^2-b^2)+10*a*b)+(4*a^4-b^4)-a^2*b^2-a*b*(a^2+b^2);
             0  0  (phi_z*L^2/3*(-10*(2*a^2-b^2)+10*a*b)+(2*(4*a^4-b^4)-2*a^2*b^2-5*L*(3*a^3-b^3)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)))/(1+phi_z)  -(phi_z*L/4*(-10*(3*a^3-b^3)-20*L*(2*a^2-b^2)+10*a*b*(a+b)+20*L*a*b)+(-10*L*(3*a^3-b^3)+3*(4*a^4-b^4)-3*a*b*(a^2+b^2)-3*a^2*b^2+10*L*a*b*(a+b)))/(1+phi_z)  5*L*(b^3+a*b^2+a^2*b-3*a^3)-((phi_z*L)^2/6*(-10*(2*a^2-b^2)+10*a*b)+phi_z/2*(2*(4*a^4-b^4)-2*a^2*b^2-5*L*(3*a^3-b^3)-2*a*b*(a^2+b^2)+5*L*a*b*(a+b)))/(1+phi_z)  phi_z*L^2/6*(-10*(2*a^2-b^2)+10*a*b)+(4*a^4-b^4)-a^2*b^2-a*b*(a^2+b^2);
             0  0  0  0  0  0;
             0  0  (phi_y*L^2/3*(-10*(a^2-2*b^2)-10*a*b)+(2*(a^4-4*b^4)+2*a^2*b^2-5*L*(a^3-3*b^3)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)))/(1+phi_y)  -(phi_y*L/4*(-10*(a^3-3*b^3)-20*L*(a^2-2*b^2)-10*a*b*(a+b)-20*L*a*b)+(-10*L*(a^3-3*b^3)+3*(a^4-4*b^4)+3*a*b*(a^2+b^2)+3*a^2*b^2-10*L*a*b*(a+b)))/(1+phi_y)  5*L*(3*b^3-a*b^2-a^2*b-a^3)-((phi_y*L)^2/6*(-10*(a^2-2*b^2)-10*a*b)+phi_y/2*(2*(a^4-4*b^4)+2*a^2*b^2-5*L*(a^3-3*b^3)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)))/(1+phi_y)  phi_y*L^2/6*(-10*(a^2-2*b^2)-10*a*b)+(a^4-4*b^4)+a^2*b^2+a*b*(a^2+b^2);
             0  0  (phi_z*L^2/3*(-10*(a^2-2*b^2)-10*a*b)+(2*(a^4-4*b^4)+2*a^2*b^2-5*L*(a^3-3*b^3)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)))/(1+phi_z)  -(phi_z*L/4*(-10*(a^3-3*b^3)-20*L*(a^2-2*b^2)-10*a*b*(a+b)-20*L*a*b)+(-10*L*(a^3-3*b^3)+3*(a^4-4*b^4)+3*a*b*(a^2+b^2)+3*a^2*b^2-10*L*a*b*(a+b)))/(1+phi_z)  5*L*(3*b^3-a*b^2-a^2*b-a^3)-((phi_z*L)^2/6*(-10*(a^2-2*b^2)-10*a*b)+phi_z/2*(2*(a^4-4*b^4)+2*a^2*b^2-5*L*(a^3-3*b^3)+2*a*b*(a^2+b^2)-5*L*a*b*(a+b)))/(1+phi_z)  phi_z*L^2/6*(-10*(a^2-2*b^2)-10*a*b)+(a^4-4*b^4)+a^2*b^2+a*b*(a^2+b^2)];
     
        for n=1:nVar
            error('Sensitivities for shear deformation in beam elements with partial DLoads are not implemented yet.')
        end
        
    end
    
    Points1 = Points(Points<=(a/L));                        % points to the left of the distributed load
    Points3 = Points(Points>=(b/L));                        % points to the right of the distributed load
    Points3 = Points3(~ismember(Points3,Points1));          
    Points2 = Points(~ismember(Points,[Points1, Points3])); % points on which the distributed load is acting
    % NeDLoad
    NeDLoad1 = zeros(length(Points1),6);
    NeDLoad2 = zeros(length(Points2),6);
    NeDLoad3 = zeros(length(Points3),6);
    for k=[2,3,4,5]
        if ~isempty(Points1), NeDLoad1(:,k) = polyval(A1(k,:),Points1); end
        if ~isempty(Points2), NeDLoad2(:,k) = polyval(A2(k,:),Points2); end
        if ~isempty(Points3), NeDLoad3(:,k) = polyval(A3(k,:),Points3); end
    end
    NeDLoad = cat(1,NeDLoad1,NeDLoad2,NeDLoad3);
       
    %dNeDLoaddx
    dNeDLoad1dx = zeros(length(Points1),6,nVar);
    dNeDLoad2dx = zeros(length(Points2),6,nVar);
    dNeDLoad3dx = zeros(length(Points3),6,nVar);
    for n=1:nVar
        for k=[2,3,4,5]
            if ~isempty(Points1), dNeDLoad1dx(:,k,n) = polyval(dA1dx(k,:,n),Points1); end
            if ~isempty(Points2), dNeDLoad2dx(:,k,n) = polyval(dA2dx(k,:,n),Points2); end
            if ~isempty(Points3), dNeDLoad3dx(:,k,n) = polyval(dA3dx(k,:,n),Points3); end
        end
    end
    dNeDLoaddx = cat(1,dNeDLoad1dx,dNeDLoad2dx,dNeDLoad3dx);
    
    % postprocessing for compatibility with the output of other 'nelcs' functions
    NeDLoad = NeDLoad/L^4;
    dNeDLoaddx = (dNeDLoaddx-NeDLoad*4*L^3.*permute(dLdx(:),[3 2 1]))/L^4;
    
end
    
