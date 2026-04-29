function [SeGCS,SeLCS,vLCS] = se_shell6(Node,Section,Material,UeGCS,Options,gcs)

%SE_SHELL6   Compute the element stresses for a shell6 element.
%
%   [SeGCS,SeLCS,vLCS] = se_shell6(Node,Section,Material,UeGCS,Options,gcs)
%   [SeGCS,SeLCS]      = se_shell6(Node,Section,Material,UeGCS,Options,gcs)
%    SeGCS             = se_shell6(Node,Section,Material,UeGCS,Options,gcs)
%   computes the element stresses in the global and the  
%   local coordinate system for the shell6 element.
%
%   Node       Node definitions           [x y z] (6 * 3)
%              Nodes should have the following order:
%              3
%              | \
%              6  5
%              |    \
%              1--4--2
%
%   Section    Section definition         [h] or [h1 h2 h3]
%              (uniform thickness or defined in corner nodes(1,2,3))
%   Material   Material definition        [E nu rho] or [Exx Eyy nuxy muxy muyz muzx theta rho]
%   UeGCS      Displacements (36 * nTimeSteps)
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
%              (in a triangular element the fourth node is zero)
%                                        [sxx syy szz sxy syz sxz]
%   SeLCS      Element stresses in LCS in corner nodes IJKL and 
%              at top/mid/bot of shell (54 * nTimeSteps)   
%                                        [sxx syy szz sxy syz sxz]
%   vLCS       Unit vectors of LCS (3 * 3)
%
%   See also ELEMSTRESS, SE_SHELL4.

% Section
h=Section(:);
if length(h)<3
    h = h(1);
elseif any(isnan(h(2:3)))
    h = h(1);
else
h = [h(1:3);(h(1)+h(2))/2;(h(2)+h(3))/2;(h(1)+h(3))/2];
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
Node=Node(1:6,1:3);
if ~ all(isfinite(Node))
    error('Not all the nodes exist.')
end

% integratiepunten
[x,H] = gaussqtri(3);

%bepaling k-factor (formule op p14-222 van ansys theorie)
A = norm(cross(Node(2,:)-Node(1,:),Node(4,:)-Node(1,:)))/2;
k=max([1.2 1+0.2*A/(25*mean(h)^2)]); 

% Constitutieve matrix
if isfield(Options,'LCSType')
LCSType = Options.LCSType;    
else
LCSType = 'element';    
end
if isfield(Options,'MatType')
MatType = Options.MatType;
else
MatType = 'isotropic';    
end    
D = cmat_shell8(MatType,Material,k);            
            

xi_eta=[0 0; 1 0; 0 1; 0.5 0; 0.5 0.5;0 0.5];
v3i = zeros(6,3);
for iNode = 1:6
    [Ni,dN_dxi,dN_deta] = sh_t6(xi_eta(iNode,1),xi_eta(iNode,2));
    Jm = [dN_dxi.'*Node(:,1) dN_dxi.'*Node(:,2) dN_dxi.'*Node(:,3);
         dN_deta.'*Node(:,1) dN_deta.'*Node(:,2) dN_deta.'*Node(:,3)];
    v3i(iNode,:) = cross(Jm(1,:),Jm(2,:));
end       

% referentie lokaal assenkruis
switch lower(LCSType)
    case 'element'
v10 = Node(2,:)-Node(1,:);
v20 = Node(3,:)-Node(1,:);
v10 = v10/norm(v10);
v30 = cross(v10,v20);
v30 = v30/norm(v30);
v20 = cross(v30,v10);
v20 = v20/norm(v20);
    case 'global'
v10 = [1 0 0];v20 = [0 1 0];v30 = [0 0 1];               
end

if nargout > 2
vLCS = [v10,v20,v30];
end

v1i = zeros(6,3);
v2i = zeros(6,3);
for iNode=1:6
    v3i(iNode,:) = v3i(iNode,:)/norm(v3i(iNode,:));
    v2i(iNode,:) = cross(v3i(iNode,:),v10);
    if norm(v2i(iNode,:))<=0.02
        v1i(iNode,:)= cross(v20,v3i(iNode,:));
        v2i(iNode,:)= cross(v3i(iNode,:),v1i(iNode,:));
    else
        v1i(iNode,:)= cross(v2i(iNode,:),v3i(iNode,:));
    end
    v1i(iNode,:) = v1i(iNode,:)/norm(v1i(iNode,:));
    v2i(iNode,:) = v2i(iNode,:)/norm(v2i(iNode,:));
end


%rotatiecomponenten transformeren naar globaal assenstelsel
T = blkdiag(eye(3), [ v2i(1,:); v1i(1,:); v3i(1,:)]);
for iNode =2:6
    T = blkdiag(T,eye(3),[ v2i(iNode,:); v1i(iNode,:); v3i(iNode,:)]);
end


SeLCSg = zeros(72,size(UeGCS,2));

for iGauss=1:size(x,1)
    xi = x(iGauss,1);
    eta = x(iGauss,2);
    
    % 2d vormfuncties
    [Ni,dN_dxi,dN_deta] = sh_t6(xi,eta);
    
    v3g = Ni.'*v3i;
    v2g = Ni.'*v2i;
    v1g = Ni.'*v1i;
     
    %transformatie matrix lokaal <-> globaal
    A = [v1g/norm(v1g); v2g/norm(v2g); v3g/norm(v3g)];
    
    % transformatie matrix uit Cook p274 en p194 Zienkiewicz deel1
    theta = [A(1:2,:).^2 A(1:2,1).*A(1:2,2) A(1:2,2).*A(1:2,3) A(1:2,3).*A(1:2,1);
             2*A(1,:).*A(2,:) A(1,1)*A(2,2)+A(1,2)*A(2,1) A(1,2)*A(2,3)+A(1,3)*A(2,2) A(1,3)*A(2,1)+A(1,1)*A(2,3);
             2*A(2,:).*A(3,:) A(2,1)*A(3,2)+A(2,2)*A(3,1) A(2,2)*A(3,3)+A(2,3)*A(3,2) A(2,3)*A(3,1)+A(2,1)*A(3,3);
             2*A(3,:).*A(1,:) A(3,1)*A(1,2)+A(3,2)*A(1,1) A(3,2)*A(1,3)+A(3,3)*A(1,2) A(3,3)*A(1,1)+A(3,1)*A(1,3)];
    
    % itereren over de hoogte (top,mid,bot)
    zeta=1; %top
    Bg = b_shell6(Ni,dN_dxi,dN_deta,zeta,Node,h,v1i,v2i,v3i,d);
    SeLCSg([1,2,4,5,6]+6*(iGauss-1),:) = D*theta*Bg*T*UeGCS;
    SeLCSg((5:6)+6*(iGauss-1),:) = 0;
    
    zeta=0; %mid
    Bg = b_shell6(Ni,dN_dxi,dN_deta,zeta,Node,h,v1i,v2i,v3i,d);
    SeLCSg([1,2,4,5,6]+6*(iGauss-1)+24,:) = D*theta*Bg*T*UeGCS;
    SeLCSg((5:6)+6*(iGauss-1)+24,:) = 1.5*SeLCSg((5:6)+6*(iGauss-1)+24,:);
    
    zeta=-1; %bot
    Bg = b_shell6(Ni,dN_dxi,dN_deta,zeta,Node,h,v1i,v2i,v3i,d);
    SeLCSg([1,2,4,5,6]+6*(iGauss-1)+48,:) = D*theta*Bg*T*UeGCS;
    SeLCSg((5:6)+6*(iGauss-1)+48,:) = 0;

end

% spanningen extrapoleren naar de knopen
extrap = inv([1-x(:,1)-x(:,2),x(:,1),x(:,2),]);
extrap = blkdiag(extrap,1);
extrap = kron(extrap,eye(6));      

Extrap = blkdiag(extrap,extrap,extrap);      
SeLCS = Extrap*SeLCSg;

SeGCS = zeros(72,size(UeGCS,2));
for iNode=1:3
    % transformatie matrix Tsigma p194 Zienkiewicz deel 1
    A = [v1i(iNode,:); v2i(iNode,:); v3i(iNode,:)]';
    theta = [A(:,1:2).^2 2*A(:,1).*A(:,2) 2*A(:,2).*A(:,3) 2*A(:,3).*A(:,1);
             A(1,1:2).*A(2,1:2) A(1,1)*A(2,2)+A(1,2)*A(2,1) A(1,2)*A(2,3)+A(1,3)*A(2,2) A(1,3)*A(2,1)+A(1,1)*A(2,3);
             A(2,1:2).*A(3,1:2) A(2,1)*A(3,2)+A(2,2)*A(3,1) A(2,2)*A(3,3)+A(2,3)*A(3,2) A(2,3)*A(3,1)+A(2,1)*A(3,3);
             A(3,1:2).*A(1,1:2) A(3,1)*A(1,2)+A(3,2)*A(1,1) A(3,2)*A(1,3)+A(3,3)*A(1,2) A(3,3)*A(1,1)+A(3,1)*A(1,3)];
   for z = 1:3; 
       SeGCS((1:6)+6*(iNode-1)+24*(z-1),:) = theta*SeLCS([1,2,4,5,6]+6*(iNode-1)+24*(z-1),:);
   end
end

if ~isempty(gcs)
    
    switch lower(gcs)
        case {'cart'}
            
        case {'cyl'}    
       a1 = [Node(1:3,1) Node(1:3,2) zeros(3,1)];
       a2 = [-Node(1:3,2) Node(1:3,1) zeros(3,1)];
       
        for iNode=1:3
        A= [a1(iNode,:)/norm(a1(iNode,:));a2(iNode,:)/norm(a2(iNode,:));0 0 1];
        theta = [A(:,:).^2 2*A(:,1).*A(:,2) 2*A(:,2).*A(:,3) 2*A(:,3).*A(:,1);
                 A(1,:).*A(2,:) A(1,1)*A(2,2)+A(1,2)*A(2,1) A(1,2)*A(2,3)+A(1,3)*A(2,2) A(1,3)*A(2,1)+A(1,1)*A(2,3);
                 A(2,:).*A(3,:) A(2,1)*A(3,2)+A(2,2)*A(3,1) A(2,2)*A(3,3)+A(2,3)*A(3,2) A(2,3)*A(3,1)+A(2,1)*A(3,3);
                 A(3,:).*A(1,:) A(3,1)*A(1,2)+A(3,2)*A(1,1) A(3,2)*A(1,3)+A(3,3)*A(1,2) A(3,3)*A(1,1)+A(3,1)*A(1,3)];
           for z = 1:3; 
           SeGCS((1:6)+6*(iNode-1)+24*(z-1),:)= theta*SeGCS((1:6)+6*(iNode-1)+24*(z-1),:);
           end
        end
   
        case {'sph'}
        a1 = Node(1:3,:);
        a2 = [-Node(1:3,2) Node(1:3,1) zeros(3,1)];
        a3 = [-Node(1:3,1).*Node(1:3,3) -Node(1:3,2).*Node(1:3,3) Node(1:3,1).^2+Node(1:3,2).^2];
        for iNode=1:3
        A= [a1(iNode,:)/norm(a1(iNode,:));a2(iNode,:)/norm(a2(iNode,:));a3(iNode,:)/norm(a3(iNode,:))];
        theta = [A(:,:).^2 2*A(:,1).*A(:,2) 2*A(:,2).*A(:,3) 2*A(:,3).*A(:,1);
                 A(1,:).*A(2,:) A(1,1)*A(2,2)+A(1,2)*A(2,1) A(1,2)*A(2,3)+A(1,3)*A(2,2) A(1,3)*A(2,1)+A(1,1)*A(2,3);
                 A(2,:).*A(3,:) A(2,1)*A(3,2)+A(2,2)*A(3,1) A(2,2)*A(3,3)+A(2,3)*A(3,2) A(2,3)*A(3,1)+A(2,1)*A(3,3);
                 A(3,:).*A(1,:) A(3,1)*A(1,2)+A(3,2)*A(1,1) A(3,2)*A(1,3)+A(3,3)*A(1,2) A(3,3)*A(1,1)+A(3,1)*A(1,3)];
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