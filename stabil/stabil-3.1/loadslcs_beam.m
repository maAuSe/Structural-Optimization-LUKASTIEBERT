function [FLCS,dFLCSdx] = loadslcs_beam(DLoadLCS,L,dDLoadLCSdx,dLdx)

%LOADSLCS_BEAM   Equivalent nodal forces for a beam element in the LCS.
%
%   FLCS = LOADSLCS_BEAM(DLoadLCS,L)
%   computes the equivalent nodal forces of a distributed load 
%   (in the local coordinate system).
%
%   [FLCS,dFLCSdx] = LOADSLCS_BEAM(DLoadLCS,L,dDLoadLCSdx,dLdx)
%   additionally computes the derivatives of the equivalent nodal forces
%   with respect to the design variables x.
%
%   DLoadLCS     Distributed loads        [n1localX; n1localY; n1localZ; ...]  (6/8 * nLC)
%   L            Beam length 
%   dDLoadLCSdx  Distributed loads derivatives	(6/8 * nLC * nVar)
%   dLdx         Beam length derivatives        (1 * nVar)
%   FLCS         Load vector                    (12 * nLC)
%   dFLCSdx      Load vector derivatives        (12 * nLC * nVar)
%
%   See also LOADS_BEAM.

% David Dooms, Wouter Dillen
% March 2008, April 2017

if nargin<3, dDLoadLCSdx = []; end
if nargin<4, dLdx = []; end

nVar = 0;
if nargout>1 && (~isempty(dDLoadLCSdx) || ~isempty(dLdx))
    nVar = max(size(dDLoadLCSdx,3),length(dLdx));
end

if nVar==0 || isempty(dDLoadLCSdx), dDLoadLCSdx = zeros(size(DLoadLCS,1),size(DLoadLCS,2),nVar); end
if nVar==0 || isempty(dLdx), dLdx = zeros(1,nVar); end

FLCS = zeros(12,size(DLoadLCS,2));
dFLCSdx = zeros(size(FLCS,1),size(FLCS,2),nVar);


% DLoad categorization
partialDLoad=-1;
if ~(nnz(DLoadLCS(1:6,:))==0)
    if size(DLoadLCS,1)==8 && ~isnan(sum(sum(DLoadLCS(7:8,:))))
        partialDLoad=1;
    else 
        partialDLoad=0;
    end
end

%% DLoad on the complete element %%
if partialDLoad==0
    
    % loads in the local x-direction
    FLCS(1,:)=L/6*(2*DLoadLCS(1,:)+DLoadLCS(4,:));
    FLCS(7,:)=L/6*(DLoadLCS(1,:)+2*DLoadLCS(4,:));

    % loads in the local y-direction
    FLCS(2,:)=L/20*(7*DLoadLCS(2,:)+3*DLoadLCS(5,:));
    FLCS(8,:)=L/20*(3*DLoadLCS(2,:)+7*DLoadLCS(5,:));
    FLCS(6,:)=L^2/60*(3*DLoadLCS(2,:)+2*DLoadLCS(5,:));
    FLCS(12,:)=-L^2/60*(2*DLoadLCS(2,:)+3*DLoadLCS(5,:));

    % loads in the local z-direction
    FLCS(3,:)=L/20*(7*DLoadLCS(3,:)+3*DLoadLCS(6,:));
    FLCS(9,:)=L/20*(3*DLoadLCS(3,:)+7*DLoadLCS(6,:));
    FLCS(5,:)=-L^2/60*(3*DLoadLCS(3,:)+2*DLoadLCS(6,:));
    FLCS(11,:)=L^2/60*(2*DLoadLCS(3,:)+3*DLoadLCS(6,:));
    
    for n=1:nVar
                    
            % load in the local x-direction derivatives
            dFLCSdx(1,:,n)=FLCS(1,:)/L*dLdx(n)+L/6*(2*dDLoadLCSdx(1,:,n)+dDLoadLCSdx(4,:,n));
            dFLCSdx(7,:,n)=FLCS(7,:)/L*dLdx(n)+L/6*(dDLoadLCSdx(1,:,n)+2*dDLoadLCSdx(4,:,n));

            % loads in the local y-direction derivatives
            dFLCSdx(2,:,n)=FLCS(2,:)/L*dLdx(n)+L/20*(7*dDLoadLCSdx(2,:,n)+3*dDLoadLCSdx(5,:,n));
            dFLCSdx(8,:,n)=FLCS(8,:)/L*dLdx(n)+L/20*(3*dDLoadLCSdx(2,:,n)+7*dDLoadLCSdx(5,:,n));
            dFLCSdx(6,:,n)=FLCS(6,:)*2/L*dLdx(n)+L^2/60*(3*dDLoadLCSdx(2,:,n)+2*dDLoadLCSdx(5,:,n));
            dFLCSdx(12,:,n)=FLCS(12,:)*2/L*dLdx(n)-L^2/60*(2*dDLoadLCSdx(2,:,n)+3*dDLoadLCSdx(5,:,n));

            % loads in the local z-direction derivatives
            dFLCSdx(3,:,n)=FLCS(3,:)/L*dLdx(n)+L/20*(7*dDLoadLCSdx(3,:,n)+3*dDLoadLCSdx(6,:,n));
            dFLCSdx(9,:,n)=FLCS(9,:)/L*dLdx(n)+L/20*(3*dDLoadLCSdx(3,:,n)+7*dDLoadLCSdx(6,:,n));
            dFLCSdx(5,:,n)=FLCS(5,:)*2/L*dLdx(n)-L^2/60*(3*dDLoadLCSdx(3,:,n)+2*dDLoadLCSdx(6,:,n));
            dFLCSdx(11,:,n)=FLCS(11,:)*2/L*dLdx(n)+L^2/60*(2*dDLoadLCSdx(3,:,n)+3*dDLoadLCSdx(6,:,n));
            
    end
    
end

%% DLoad on part of the element %%
if partialDLoad==1
    
    % j = DLoadLCS column index containing the DLoad information (should be only one)
    j = find(sum(abs(DLoadLCS(1:6,:)),1));
    if length(j) > 1
        error('Wrong DLoads structure in dimension 3; make sure to use MULTDLOADS when combining multiple load cases.');
    end
    
    % set load starting & ending point
    a = DLoadLCS(7,j);
    b = DLoadLCS(8,j);

    % loads in the local x-direction
    q_a = DLoadLCS(1,j);
    q_b = DLoadLCS(4,j);
    FLCS(1,j) = 1/(b-a)*((a^2-2*a*b+b^2)/2*(q_b+q_a)) - 1/(b-a)/L*((a^3-3*a*b^2+2*b^3)/6*q_b+(2*a^3-3*a^2*b+b^3)/6*q_a);    % Q_ix
    FLCS(7,j) = 1/(b-a)/L*((a^3-3*a*b^2+2*b^3)/6*q_b+(2*a^3-3*a^2*b+b^3)/6*q_a);                                            % Q_jx
    for n=1:nVar
        dadx = dDLoadLCSdx(7,j,n);
        dbdx = dDLoadLCSdx(8,j,n);
        dqadx = dDLoadLCSdx(1,j,n);
        dqbdx = dDLoadLCSdx(4,j,n);
        dFLCSdx(1,j,n) = (((dadx-dbdx)/(b-a)^2)*((a^2-2*a*b+b^2)/2*(q_b+q_a)) + 1/(b-a)*((a-b)*(dadx-dbdx)*(q_b+q_a)+(a^2-2*a*b+b^2)/2*(dqbdx+dqadx))) - (((dadx-dbdx)/(b-a)^2/L-1/(b-a)*dLdx(n)/L^2)*((a^3-3*a*b^2+2*b^3)/6*q_b+(2*a^3-3*a^2*b+b^3)/6*q_a) + 1/(b-a)/L*(((a^2-b^2)/2*dadx-b*(a-b)*dbdx)*q_b+(a^3-3*a*b^2+2*b^3)/6*dqbdx+(a*(a-b)*dadx-(a^2-b^2)/2*dbdx)*q_a+(2*a^3-3*a^2*b+b^3)/6*dqadx));   % dQ_ixdx
        dFLCSdx(7,j,n) = (((dadx-dbdx)/(b-a)^2/L-1/(b-a)*dLdx(n)/L^2)*((a^3-3*a*b^2+2*b^3)/6*q_b+(2*a^3-3*a^2*b+b^3)/6*q_a) + 1/(b-a)/L*(((a^2-b^2)/2*dadx-b*(a-b)*dbdx)*q_b+(a^3-3*a*b^2+2*b^3)/6*dqbdx+(a*(a-b)*dadx-(a^2-b^2)/2*dbdx)*q_a+(2*a^3-3*a^2*b+b^3)/6*dqadx));                                                                                                                                   % dQ_jxdx
    end

    % loads in the local y-direction
    q_a = DLoadLCS(2,j);
    q_b = DLoadLCS(5,j);
    FLCS(2,j) = 1/(b-a)*((a^2-2*a*b+b^2)/2*(q_b+q_a)) - 3/(b-a)/L^2*((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a) + 2/(b-a)/L^3*((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a);              % Q_iy
    FLCS(8,j) = 3/(b-a)/L^2*((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a) - 2/(b-a)/L^3*((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a);                                                      % Q_jy
    FLCS(6,j) = 1/(b-a)*(((a^3-3*a*b^2+2*b^3)/6*q_b+(2*a^3-3*a^2*b+b^3)/6*q_a) - ((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a)*2/L + ((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a)/L^2);    % M_iz
    FLCS(12,j) = - 1/(b-a)*(((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a)/L - ((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a)/L^2);                                                           % M_jz
    for n=1:nVar
        dadx = dDLoadLCSdx(7,j,n);
        dbdx = dDLoadLCSdx(8,j,n);
        dqadx = dDLoadLCSdx(2,j,n);
        dqbdx = dDLoadLCSdx(5,j,n);
        dFLCSdx(2,j,n) = (((dadx-dbdx)/(b-a)^2)*((a^2-2*a*b+b^2)/2*(q_b+q_a)) + 1/(b-a)*((a-b)*(dadx-dbdx)*(q_b+q_a)+(a^2-2*a*b+b^2)/2*(dqbdx+dqadx))) - (((dadx-dbdx)/(b-a)^2*3/L^2-1/(b-a)*6*dLdx(n)/L^3)*((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a) + 3/(b-a)/L^2*(((a^3-b^3)/3*dadx-b^2*(a-b)*dbdx)*q_b+(a^4-4*a*b^3+3*b^4)/12*dqbdx+(a^2*(a-b)*dadx-(a^3-b^3)/3*dbdx)*q_a+(3*a^4-4*a^3*b+b^4)/12*dqadx)) + (((dadx-dbdx)/(b-a)^2*2/L^3-1/(b-a)*6*dLdx(n)/L^4)*((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a) + 2/(b-a)/L^3*(((a^4-b^4)/4*dadx-b^3*(a-b)*dbdx)*q_b+(a^5-5*a*b^4+4*b^5)/20*dqbdx+(a^3*(a-b)*dadx-(a^4-b^4)/4*dbdx)*q_a+(4*a^5-5*a^4*b+b^5)/20*dqadx));     % dQ_iydx
        dFLCSdx(8,j,n) = (((dadx-dbdx)/(b-a)^2*3/L^2-1/(b-a)*6*dLdx(n)/L^3)*((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a) + 3/(b-a)/L^2*(((a^3-b^3)/3*dadx-b^2*(a-b)*dbdx)*q_b+(a^4-4*a*b^3+3*b^4)/12*dqbdx+(a^2*(a-b)*dadx-(a^3-b^3)/3*dbdx)*q_a+(3*a^4-4*a^3*b+b^4)/12*dqadx)) - (((dadx-dbdx)/(b-a)^2*2/L^3-1/(b-a)*6*dLdx(n)/L^4)*((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a) + 2/(b-a)/L^3*(((a^4-b^4)/4*dadx-b^3*(a-b)*dbdx)*q_b+(a^5-5*a*b^4+4*b^5)/20*dqbdx+(a^3*(a-b)*dadx-(a^4-b^4)/4*dbdx)*q_a+(4*a^5-5*a^4*b+b^5)/20*dqadx));                                                                                                                                     % dQ_jydx
        dFLCSdx(6,j,n) = ((dadx-dbdx)/(b-a)^2)*(FLCS(6,j)*(b-a)) + 1/(b-a)*((((a^2-b^2)/2*dadx-b*(a-b)*dbdx)*q_b+(a^3-3*a*b^2+2*b^3)/6*dqbdx+(a*(a-b)*dadx-(a^2-b^2)/2*dbdx)*q_a+(2*a^3-3*a^2*b+b^3)/6*dqadx) - (((a^3-b^3)/3*dadx-b^2*(a-b)*dbdx)*q_b+(a^4-4*a*b^3+3*b^4)/12*dqbdx+(a^2*(a-b)*dadx-(a^3-b^3)/3*dbdx)*q_a+(3*a^4-4*a^3*b+b^4)/12*dqadx)*2/L + ((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a)*2*dLdx(n)/L^2 + (((a^4-b^4)/4*dadx-b^3*(a-b)*dbdx)*q_b+(a^5-5*a*b^4+4*b^5)/20*dqbdx+(a^3*(a-b)*dadx-(a^4-b^4)/4*dbdx)*q_a+(4*a^5-5*a^4*b+b^5)/20*dqadx)/L^2 - ((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a)*2*dLdx(n)/L^3);                                         % dM_izdx
        dFLCSdx(12,j,n) = - ((dadx-dbdx)/(b-a)^2)*(-FLCS(12,j)*(b-a)) - 1/(b-a)*((((a^3-b^3)/3*dadx-b^2*(a-b)*dbdx)*q_b+(a^4-4*a*b^3+3*b^4)/12*dqbdx+(a^2*(a-b)*dadx-(a^3-b^3)/3*dbdx)*q_a+(3*a^4-4*a^3*b+b^4)/12*dqadx)/L - ((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a)*dLdx(n)/L^2 - (((a^4-b^4)/4*dadx-b^3*(a-b)*dbdx)*q_b+(a^5-5*a*b^4+4*b^5)/20*dqbdx+(a^3*(a-b)*dadx-(a^4-b^4)/4*dbdx)*q_a+(4*a^5-5*a^4*b+b^5)/20*dqadx)/L^2 + ((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a)*2*dLdx(n)/L^3);                                                                                                                                                                            % dM_jzdx
    end
    
    % loads in the local z-direction
    q_a = DLoadLCS(3,j);
    q_b = DLoadLCS(6,j);
    FLCS(3,j) = 1/(b-a)*((a^2-2*a*b+b^2)/2*(q_b+q_a)) - 3/(b-a)/L^2*((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a) + 2/(b-a)/L^3*((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a);              % Q_iz
    FLCS(9,j) = 3/(b-a)/L^2*((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a) - 2/(b-a)/L^3*((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a);                                                      % Q_jz
    FLCS(5,j) = - 1/(b-a)*(((a^3-3*a*b^2+2*b^3)/6*q_b+(2*a^3-3*a^2*b+b^3)/6*q_a) - ((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a)*2/L + ((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a)/L^2);  % M_iy
    FLCS(11,j) = 1/(b-a)*(((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a)/L - ((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a)/L^2);                                                             % M_jz
    for n=1:nVar
        dadx = dDLoadLCSdx(7,j,n);
        dbdx = dDLoadLCSdx(8,j,n);
        dqadx = dDLoadLCSdx(3,j,n);
        dqbdx = dDLoadLCSdx(6,j,n);
        dFLCSdx(3,j,n) = (((dadx-dbdx)/(b-a)^2)*((a^2-2*a*b+b^2)/2*(q_b+q_a)) + 1/(b-a)*((a-b)*(dadx-dbdx)*(q_b+q_a)+(a^2-2*a*b+b^2)/2*(dqbdx+dqadx))) - (((dadx-dbdx)/(b-a)^2*3/L^2-1/(b-a)*6*dLdx(n)/L^3)*((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a) + 3/(b-a)/L^2*(((a^3-b^3)/3*dadx-b^2*(a-b)*dbdx)*q_b+(a^4-4*a*b^3+3*b^4)/12*dqbdx+(a^2*(a-b)*dadx-(a^3-b^3)/3*dbdx)*q_a+(3*a^4-4*a^3*b+b^4)/12*dqadx)) + (((dadx-dbdx)/(b-a)^2*2/L^3-1/(b-a)*6*dLdx(n)/L^4)*((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a) + 2/(b-a)/L^3*(((a^4-b^4)/4*dadx-b^3*(a-b)*dbdx)*q_b+(a^5-5*a*b^4+4*b^5)/20*dqbdx+(a^3*(a-b)*dadx-(a^4-b^4)/4*dbdx)*q_a+(4*a^5-5*a^4*b+b^5)/20*dqadx));     % dQ_izdx
        dFLCSdx(9,j,n) = (((dadx-dbdx)/(b-a)^2*3/L^2-1/(b-a)*6*dLdx(n)/L^3)*((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a) + 3/(b-a)/L^2*(((a^3-b^3)/3*dadx-b^2*(a-b)*dbdx)*q_b+(a^4-4*a*b^3+3*b^4)/12*dqbdx+(a^2*(a-b)*dadx-(a^3-b^3)/3*dbdx)*q_a+(3*a^4-4*a^3*b+b^4)/12*dqadx)) - (((dadx-dbdx)/(b-a)^2*2/L^3-1/(b-a)*6*dLdx(n)/L^4)*((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a) + 2/(b-a)/L^3*(((a^4-b^4)/4*dadx-b^3*(a-b)*dbdx)*q_b+(a^5-5*a*b^4+4*b^5)/20*dqbdx+(a^3*(a-b)*dadx-(a^4-b^4)/4*dbdx)*q_a+(4*a^5-5*a^4*b+b^5)/20*dqadx));                                                                                                                                     % dQ_jzdx
        dFLCSdx(5,j,n) = - ((dadx-dbdx)/(b-a)^2)*(-FLCS(5,j)*(b-a)) - 1/(b-a)*((((a^2-b^2)/2*dadx-b*(a-b)*dbdx)*q_b+(a^3-3*a*b^2+2*b^3)/6*dqbdx+(a*(a-b)*dadx-(a^2-b^2)/2*dbdx)*q_a+(2*a^3-3*a^2*b+b^3)/6*dqadx) - (((a^3-b^3)/3*dadx-b^2*(a-b)*dbdx)*q_b+(a^4-4*a*b^3+3*b^4)/12*dqbdx+(a^2*(a-b)*dadx-(a^3-b^3)/3*dbdx)*q_a+(3*a^4-4*a^3*b+b^4)/12*dqadx)*2/L + ((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a)*2*dLdx(n)/L^2 + (((a^4-b^4)/4*dadx-b^3*(a-b)*dbdx)*q_b+(a^5-5*a*b^4+4*b^5)/20*dqbdx+(a^3*(a-b)*dadx-(a^4-b^4)/4*dbdx)*q_a+(4*a^5-5*a^4*b+b^5)/20*dqadx)/L^2 - ((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a)*2*dLdx(n)/L^3);                                      % dM_iydx
        dFLCSdx(11,j,n) = ((dadx-dbdx)/(b-a)^2)*(FLCS(11,j)*(b-a)) + 1/(b-a)*((((a^3-b^3)/3*dadx-b^2*(a-b)*dbdx)*q_b+(a^4-4*a*b^3+3*b^4)/12*dqbdx+(a^2*(a-b)*dadx-(a^3-b^3)/3*dbdx)*q_a+(3*a^4-4*a^3*b+b^4)/12*dqadx)/L - ((a^4-4*a*b^3+3*b^4)/12*q_b+(3*a^4-4*a^3*b+b^4)/12*q_a)*dLdx(n)/L^2 - (((a^4-b^4)/4*dadx-b^3*(a-b)*dbdx)*q_b+(a^5-5*a*b^4+4*b^5)/20*dqbdx+(a^3*(a-b)*dadx-(a^4-b^4)/4*dbdx)*q_a+(4*a^5-5*a^4*b+b^5)/20*dqadx)/L^2 + ((a^5-5*a*b^4+4*b^5)/20*q_b+(4*a^5-5*a^4*b+b^5)/20*q_a)*2*dLdx(n)/L^3);                                                                                                                                                                               % dM_jydx
    end
end
