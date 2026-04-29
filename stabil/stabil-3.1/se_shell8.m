function [SeGCS,SeLCS,vLCS] = se_shell8(Node,Section,Material,UeGCS,Options,gcs)

%SE_SHELL8   Compute the element stresses for a shell8 element.
%
%   [SeGCS,SeLCS,vLCS] = se_shell8(Node,Section,Material,UeGCS,Options,gcs)
%   [SeGCS,SeLCS]      = se_shell8(Node,Section,Material,UeGCS,Options,gcs)
%    SeGCS             = se_shell8(Node,Section,Material,UeGCS,Options,gcs)
%   computes the element stresses in the global and the
%   local coordinate system for the shell8 element.
%
%   Node       Node definitions           [x y z] (8 * 3)
%              Nodes should have the following order:
%              4----7----3
%              |         |
%              8         6
%              |         |
%              1----5----2
%
%   Section    Section definition         [h] or [h1 h2 h3 h4]
%              (uniform thickness or defined in corner nodes(1,2,3,4))
%   Material   Material definition        [E nu rho] or [Exx Eyy nuxy muxy muyz muzx theta rho]
%   UeGCS      Displacements (48 * nTimeSteps)
%   Options    Element options struct. Fields:
%              -LCSType: determine the reference local element
%                        coordinate system. Values:
%                        'element' (default) or 'global'
%              -MatType: 'isotropic' (default) or 'orthotropic'
%              -Offset: nodal offset from shell midplane. Values:
%                      'top', 'mid' (default), 'bot' or numerical value
%   GCS        Global coordinate system in which stresses are returned
%              'cart'|'cyl'|'sph'
%   SeGCS      Element stresses in GCS in corner nodes IJKL and
%              at top/mid/bot of shell (72 * nTimeSteps)
%              72 = 6 stress comp. * 4 nodes * 3 locations
%                                        [sxx syy szz sxy syz sxz]
%   SeLCS      Element stresses in LCS in corner nodes IJKL and
%              at top/mid/bot of shell (72 * nTimeSteps)
%                                        [sxx syy szz sxy syz sxz]
%   vLCS       Unit vectors of LCS (3 * 3)
%
%   See also ELEMSTRESS, SE_SHELL4.

% Miche Jansen
% 2010

% Section
h=Section(:);
if length(h)<4
    h = h(1);
elseif any(isnan(h(2:4)))
    h = h(1);
else
    h = [h(1:4);(h(1)+h(2))/2;(h(2)+h(3))/2;(h(3)+h(4))/2;(h(1)+h(4))/2];
end
% offset
if isfield(Options,'Offset')
    if isnumeric(Options.Offset)
        d = Options.Offset;
    else
        switch lower(Options.Offset)
            case 'top'; d = h/2;
            case 'mid'; d = 0;
            case 'bot'; d = -h/2;
        end
    end
else
    d = 0;
end

% Check nodes
Node=Node(1:8,1:3);
if ~ all(isfinite(Node))
    error('Not all the nodes exist.')
end

% integratiepunten
[x,H] = gaussq(2);


A = norm(cross(Node(2,:)-Node(1,:),Node(4,:)-Node(1,:)));
k=max([1.2 1+0.2*A/(25*mean(h)^2)]);

% Constitutieve matrix
if isfield(Options,'MatType')
    MatType = Options.MatType;
else
    MatType = 'isotropic';
end
D = cmat_shell8(MatType,Material,k);


xi_eta = [ -1 -1;1  -1;1 1; -1 1;0 -1; 1 0; 0 1; -1 0];
v3i = zeros(8,3);
for iNode = 1:8
    [Ni,dN_dxi,dN_deta] = sh_qs8(xi_eta(iNode,1),xi_eta(iNode,2));
    Jm = [dN_dxi.'*Node(:,1) dN_dxi.'*Node(:,2) dN_dxi.'*Node(:,3);
        dN_deta.'*Node(:,1) dN_deta.'*Node(:,2) dN_deta.'*Node(:,3)];
    v3i(iNode,:) = cross(Jm(1,:),Jm(2,:));
end

% referentie lokaal assenkruis
v0 = trans_shell8(Node,Options);
if nargout > 2
    vLCS = v0.'; vLCS = vLCS(:).';
end

% lokaal assenkruis per knoop --> projectie-methode COOK p579
v1i = zeros(8,3);
v2i = zeros(8,3);
for iNode=1:8
    v3i(iNode,:) = v3i(iNode,:)/norm(v3i(iNode,:));
    v2i(iNode,:) = cross(v3i(iNode,:),v0(1,:));
    if norm(v2i(iNode,:))<=0.02
        v1i(iNode,:)= cross(v0(2,:),v3i(iNode,:));
        v2i(iNode,:)= cross(v3i(iNode,:),v1i(iNode,:));
    else
        v1i(iNode,:)= cross(v2i(iNode,:),v3i(iNode,:));
    end
    v1i(iNode,:) = v1i(iNode,:)/norm(v1i(iNode,:));
    v2i(iNode,:) = v2i(iNode,:)/norm(v2i(iNode,:));
end


%rotatiecomponenten transformeren naar globaal assenstelsel
T = blkdiag(eye(3), [ v2i(1,:); v1i(1,:); v3i(1,:)]);
for iNode =2:8
    T = blkdiag(T,eye(3),[ v2i(iNode,:); v1i(iNode,:); v3i(iNode,:)]);
end

UeGLCS = T*UeGCS; % displacement vector with translations in GCS and rotations in LCS
SeLCSg = zeros(72,size(UeGCS,2));

for iGauss=1:size(x,1)
    xi = x(iGauss,1);
    eta = x(iGauss,2);
    
    % 2d vormfuncties
    [Ni,dN_dxi,dN_deta] = sh_qs8(xi,eta);
    
    %bepaling lokale z-richting bij het gausspunt
    %     Jm = [dN_dxi.'*Node(:,1) dN_dxi.'*Node(:,2) dN_dxi.'*Node(:,3);
    %          dN_deta.'*Node(:,1) dN_deta.'*Node(:,2) dN_deta.'*Node(:,3)];
    %     v3g=cross(Jm(1,:),Jm(2,:));
    %     [mini,Ig]=min(v3g);
    %     v1g = cross(Ei(Ig,:),v3g);
    %     v2g = cross(v3g,v1g);
    v3g = Ni.'*v3i;
    v2g = Ni.'*v2i;
    v1g = Ni.'*v1i;
    
    %transformatie matrix lokaal <-> globaal
    A = [v1g/norm(v1g); v2g/norm(v2g); v3g/norm(v3g)];
    
    % transformatie matrix uit Cook p274 en p194 Zienkiewicz deel1
    theta = vtrans_solid(A,'strain');
    % constitutive matrix for epsilon_gcs -> sigma_lcs
    D2 = D*theta([1,2,4,5,6],:);
    
    % itereren over de hoogte (top,mid,bot)
    zeta=1; %top
    Bg = b_shell8(Ni,dN_dxi,dN_deta,zeta,Node,h,v1i,v2i,v3i,d);
    SeLCSg([1,2,4,5,6]+6*(iGauss-1),:) = D2*Bg*UeGLCS;
    SeLCSg((5:6)+6*(iGauss-1),:) = 0;
    %SeLCSg(4:5,:) = 0; % schuifspanningen =0 (geen sz; daarom 4:5)
    %SeGCSg((1:6)+6*(iGauss-1),:) = theta.'*SeLCSg;
    
    zeta=0; %mid
    Bg = b_shell8(Ni,dN_dxi,dN_deta,zeta,Node,h,v1i,v2i,v3i,d);
    SeLCSg([1,2,4,5,6]+6*(iGauss-1)+24,:) = D2*Bg*UeGLCS;
    SeLCSg((5:6)+6*(iGauss-1)+24,:) = 1.5*SeLCSg((5:6)+6*(iGauss-1)+24,:);
    %SeLCSg(4:5,:) = 1.5*SeLCSg(4:5,:);
    %SeGCSg((1:6)+6*(iGauss-1)+24,:) = theta.'*SeLCSg;
    
    zeta=-1; %bot
    Bg = b_shell8(Ni,dN_dxi,dN_deta,zeta,Node,h,v1i,v2i,v3i,d);
    SeLCSg([1,2,4,5,6]+6*(iGauss-1)+48,:) = D2*Bg*UeGLCS;
    SeLCSg((5:6)+6*(iGauss-1)+48,:) = 0;
    %SeLCSg(4:5,:) = 0; % schuifspanningen =0
    %SeGCSg((1:6)+6*(iGauss-1)+48,:) = theta.'*SeLCSg;
    
end

% spanningen extrapoleren naar de knopen
g=1/abs(x(1,1))/2;
extrap = [(1+g)*eye(6) -0.5*eye(6) (1-g)*eye(6) -0.5*eye(6);...
    -0.5*eye(6) (1+g)*eye(6) -0.5*eye(6) (1-g)*eye(6);...
    (1-g)*eye(6) -0.5*eye(6) (1+g)*eye(6) -0.5*eye(6);...
    -0.5*eye(6) (1-g)*eye(6) -0.5*eye(6) (1+g)*eye(6)];
Extrap = blkdiag(extrap,extrap,extrap);
SeLCS = Extrap*SeLCSg;

SeGCS = zeros(72,size(UeGCS,2));
for iNode=1:4
    % transformatie matrix Tsigma p194 Zienkiewicz deel 1
    A = [v1i(iNode,:); v2i(iNode,:); v3i(iNode,:)]'; % transformation matrix LCS -> GCS
    theta = vtrans_solid(A);
    for z = 1:3;
        SeGCS((1:6)+6*(iNode-1)+24*(z-1),:) = theta*SeLCS((1:6)+6*(iNode-1)+24*(z-1),:);
        %SeLCS([1,2,4,5,6]+6*(iNode-1)+24*(z-1),:)= theta*SeGCS((1:6)+6*(iNode-1)+24*(z-1),:);
    end
end

if ~isempty(gcs)
    
    switch lower(gcs)
        case {'cart'}
            
        case {'cyl'}
            a1 = [Node(1:4,1) Node(1:4,2) zeros(4,1)];
            a2 = [-Node(1:4,2) Node(1:4,1) zeros(4,1)];
            
            for iNode=1:4
                A= [a1(iNode,:)/norm(a1(iNode,:));a2(iNode,:)/norm(a2(iNode,:));0 0 1];
                theta = vtrans_solid(A);
                for z = 1:3;
                    SeGCS((1:6)+6*(iNode-1)+24*(z-1),:)= theta*SeGCS((1:6)+6*(iNode-1)+24*(z-1),:);
                end
            end
            
        case {'sph'}
            a1 = Node(1:4,:);
            a2 = [-Node(1:4,2) Node(1:4,1) zeros(4,1)];
            a3 = [-Node(1:4,1).*Node(1:4,3) -Node(1:4,2).*Node(1:4,3) Node(1:4,1).^2+Node(1:4,2).^2];
            for iNode=1:4
                A= [a1(iNode,:)/norm(a1(iNode,:));a2(iNode,:)/norm(a2(iNode,:));a3(iNode,:)/norm(a3(iNode,:))];
                theta = vtrans_solid(A);
                for z = 1:3;
                    SeGCS((1:6)+6*(iNode-1)+24*(z-1),:)= theta*SeGCS((1:6)+6*(iNode-1)+24*(z-1),:);
                end
            end
            
            
    end
end

end

function z = cross(x,y)
z = x([2 1 1]).*y([3 3 2])-y([2 1 1]).*x([3 3 2]);
z(2) = -z(2);
end
