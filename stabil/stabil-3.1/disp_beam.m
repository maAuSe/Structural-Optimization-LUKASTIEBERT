function [Ax,Ay,Az,B,Cx,Cy,Cz,dAxdx,dAydx,dAzdx,dCxdx,dCydx,dCzdx] = disp_beam(Nodes,Elements,DOF,DLoads,Sections,Materials,Points,dNodesdx,dDLoadsdx,dSectionsdx)

%DISP_BEAM   Return matrices to compute the displacements of the deformed beams.
%
%   [Ax,Ay,Az,B,Cx,Cy,Cz] = DISP_BEAM(Nodes,Elements,DOF,DLoads,Sections,Materials,Points)       
%   [Ax,Ay,Az,B,Cx,Cy,Cz] = DISP_BEAM(Nodes,Elements,DOF,DLoads,Sections,Materials)
%   [Ax,Ay,Az,B] = DISP_BEAM(Nodes,Elements,DOF,[],Sections,Materials) 
%   [Ax,Ay,Az,B] = DISP_BEAM(Nodes,Elements,DOF)
%       returns the matrices to compute the displacements of the deformed
%       beams. The coordinates of the specified points along the deformed
%       beam elements are computed using X=Ax*U+Cx*DLoad+B(:,1); 
%       Y=Ay*U+Cy*DLoad+B(:,2) and Z=Az*U+Cz*DLoad+B(:,3). The matrices 
%       Cx,Cy and Cz superimpose the displacements that occur due to the 
%       distributed loads if all nodes are fixed. 
%
%   [Ax,Ay,Az,B,Cx,Cy,Cz,dAxdx,dAydx,dAzdx,dCxdx,dCydx,dCzdx] 
%           = DISP_BEAM(Nodes,Elements,DOF,DLoads,Sections,Materials,Points,dNodesdx,
%                                                                 dDLoadsdx(,dSectionsdx)) 
%       additionally computes the derivatives of the displacements with
%       respect to the design variables x.
%
%   Nodes      Node definitions         [NodID x y z]
%   Elements   Element definitions      [EltID TypID SecID MatID n1 n2 ...]
%   DOF        Degrees of freedom       (nDOF * 1)
%   DLoads     Distributed loads        [EltID n1globalX n1globalY n1globalZ ...]
%                                 (use an empty array [] when shear deformation is 
%                                                   considered but no DLoads are present) 
%   Sections   Section definitions   	[SecID SecProp1 SecProp2 ...]
%   Materials  Material definitions    	[MatID MatProp1 MatProp2 ... ]
%   Points     Points in the local coordinate system  (1 * nPoints)
%   dNodesdx   Node definitions derivatives   (SIZE(Node) * nVar)
%   dDLoadsdx  Distributed loads derivatives  (SIZE(DLoad) * nVar)
%   Ax         Matrix to compute the x-coordinates of the deformations 
%   Ay         Matrix to compute the y-coordinates of the deformations 
%   Az         Matrix to compute the z-coordinates of the deformations 
%   B          Matrix which contains the x-, y- and z-coordinates of the 
%              undeformed structure
%   Cx         Matrix to compute the x-coordinates of the deformations 
%   Cy         Matrix to compute the y-coordinates of the deformations 
%   Cz         Matrix to compute the z-coordinates of the deformations 
%   dAxdx, dAydx, dAzdx, dCxdx, dCydx, dCzdx    
%       Derivatives of the matrices to compute the coordinates of the
%       interpolation points after deformation
%
%   See also DISP_TRUSS, PLOTDISP, NELCS_BEAM, NEDLOADLCS_BEAM.

% David Dooms, Wouter Dillen
% September 2008, July 2017

if nargin<4, DLoads = []; end
if nargin<5, Sections = []; end
if nargin<6, Materials = []; end
if nargin<7, Points = []; end
if nargin<8, dNodesdx = []; end
if nargin<9, dDLoadsdx = []; end
if nargin<10, dSectionsdx = []; end

nVar = 0;
if nargout>7 && (~isempty(dNodesdx) || ~isempty(dDLoadsdx) || ~isempty(dSectionsdx))
    nVar = max([size(dNodesdx,3),size(dDLoadsdx,4),size(dSectionsdx,3)]);
end

if isempty(Points)
    nPoints=21;
    Points = linspace(0,1,nPoints);
else
    Points = Points(:);
    nPoints=length(Points);
end
if nVar==0 || isempty(dNodesdx), dNodesdx = zeros([size(Nodes),nVar]); end
if nVar==0 || isempty(dDLoadsdx), dDLoadsdx = zeros(size(DLoads,1),size(DLoads,2),size(DLoads,3),nVar); end
if nVar==0 || isempty(dSectionsdx), dSectionsdx = zeros([size(Sections),nVar]); end


nElem=size(Elements,1);
DOF=DOF(:);
nDOF=length(DOF);

Ax=sparse([],[],[],nElem*(nPoints+1),nDOF,nElem*nPoints*12+nElem);
Ay=sparse([],[],[],nElem*(nPoints+1),nDOF,nElem*nPoints*12+nElem);
Az=sparse([],[],[],nElem*(nPoints+1),nDOF,nElem*nPoints*12+nElem);
B=zeros(nElem*(nPoints+1),3);
%
dAxdx=zeros(nElem*(nPoints+1),nDOF,nVar);
dAydx=zeros(nElem*(nPoints+1),nDOF,nVar);
dAzdx=zeros(nElem*(nPoints+1),nDOF,nVar);
if ~isempty(DLoads)  % en als DLoad leeg is?  materiaal en sectie ontbreken? 
    EltIDDLoad=DLoads(:,1);
    nDLoad=size(EltIDDLoad,1);
    Cx=sparse([],[],[],nElem*(nPoints+1),nDLoad*6,nDLoad*nPoints*6);
    Cy=sparse([],[],[],nElem*(nPoints+1),nDLoad*6,nDLoad*nPoints*6);
    Cz=sparse([],[],[],nElem*(nPoints+1),nDLoad*6,nDLoad*nPoints*6);
    %
    dCxdx=zeros(nElem*(nPoints+1),nDLoad*6,nVar);
    dCydx=zeros(nElem*(nPoints+1),nDLoad*6,nVar);
    dCzdx=zeros(nElem*(nPoints+1),nDLoad*6,nVar);
elseif nargout>4
    Cx = [];
    Cy = [];
    Cz = [];
end
    

% When there is no shear deformation
if isempty(Sections) && isempty(Materials)
    % if only one of them is empty, the code will throw an error later on
    NeLCS=nelcs_beam(Points);
end

% Loop over the elements
for iElem=1:nElem
    
    % Nodes
    NodeNum=Elements(iElem,5:end);
    Node=zeros(length(NodeNum),3);
    dNodedx=zeros(length(NodeNum),3,nVar);
    
    for iNode=1:length(NodeNum)
        loc=find(Nodes(:,1)==NodeNum(1,iNode));
        if isempty(loc)
            Node(iNode,:)=NaN;
            dNodedx(iNode,:,:)=NaN;
        elseif length(loc)>1
            error('Node %i is multiply defined.',NodeNum(1,iNode))
        else
            Node(iNode,:)=Nodes(loc,2:end);
            dNodedx(iNode,:,:)=dNodesdx(loc,2:end,:);
        end
    end
    
    dofelem=dof_beam(NodeNum);
    C=selectdof(DOF,dofelem);
    
    % Transform displacements from global to local coordinate system
    [t,dtdx]=trans_beam(Node,dNodedx);
    
    % Element length
    L=norm(Node(2,:)-Node(1,:)); 
    dLdx=(permute(dNodedx(2,:,:)-dNodedx(1,:,:),[3 2 1])*(Node(2,:)-Node(1,:)).')/L;
    
    % When there is DLoads and/or shear deformation
    if ~isempty(DLoads) || ~isempty(Sections) || ~isempty(Materials)
        SecID=Elements(iElem,3);
        loc=find(Sections(:,1)==SecID);
        if isempty(loc)
            error('Section %i is not defined.',SecID)
        elseif length(loc)>1
            error('Section %i is multiply defined.',SecID)
        else
            A=Sections(loc,2);
            ky=Sections(loc,3);
            kz=Sections(loc,4);
            Iyy=Sections(loc,6);
            Izz=Sections(loc,7);
            dAdx=dSectionsdx(loc,2,:);
            dkydx=dSectionsdx(loc,3,:);
            dkzdx=dSectionsdx(loc,4,:);
            dIyydx=dSectionsdx(loc,6,:);
            dIzzdx=dSectionsdx(loc,7,:);
        end
  
        MatID=Elements(iElem,4);
        loc=find(Materials(:,1)==MatID);
        if isempty(loc)
            error('Material %i is not defined.',MatID)
        elseif length(loc)>1
            error('Material %i is multiply defined.',MatID)
        else
            E=Materials(loc,2);
            nu=Materials(loc,3);
        end

        phi_y=24*(1+nu)*Izz/(A*(ky+1e-250)*L^2);
        phi_z=24*(1+nu)*Iyy/(A*(kz+1e-250)*L^2);
        NeLCS = nelcs_beam(Points,phi_y,phi_z);   
        
        if nVar>0
            if isfinite(ky)  &&  any([dAdx(:); dIzzdx(:); dLdx(:); dkydx(:)]) > 0
                % dphi_ydx is not equal to zero
                error('Sensitivities for Timoshenko beam elements have not been implemented yet.')
            end
            if isfinite(kz)  &&  any([dAdx(:); dIyydx(:); dLdx(:); dkzdx(:)]) > 0
                % dphi_zdx is not equal to zero
                error('Sensitivities for Timoshenko beam elements have not been implemented yet.')
            end
        end
     
        % When there is DLoads [calculate Cx Cy Cz]
        if ~isempty(DLoads)  
            
            EltID=Elements(iElem,1);
            iDLoad=find(EltIDDLoad==EltID);
            for i=1:length(iDLoad) % loop over all DLoads on the element
                
                irow = iDLoad(i);
                if size(DLoads,2)==9 && ~isnan(sum(sum(DLoads(irow,8:9,:),3)))
                    % partial DLoad
                    j = find(sum(abs(DLoads(irow,2:7,:)),2)~=0);
                    if isempty(j) % trivial partial DLoad of 0
                        a = 0;
                        b = L;
                        dadx = zeros(size(dLdx));
                        dbdx = zeros(size(dLdx));
                    else
                        if length(j) > 1, error('Wrong DLoads structure in dimension 3; make sure to use MULTDLOADS when combining multiple load cases.'); end
                        a = DLoads(irow,8,j);
                        b = DLoads(irow,9,j);
                        dadx = permute(dDLoadsdx(irow,8,j,:),[4 1 2 3]);
                        dbdx = permute(dDLoadsdx(irow,9,j,:),[4 1 2 3]);
                    end
                    [f1,df1dx] = nedloadlcs_beam(Points,phi_y,phi_z,a,b,L,dadx,dbdx,dLdx);
                else
                    % complete DLoad
                    f1 = nedloadlcs_beam(Points,phi_y,phi_z);
                    df1dx = zeros([size(f1),nVar]);
                end
                NeDLoad = f1*L^4;
                dNeDLoaddx = df1dx*L^4 + f1*4*L^3.*permute(dLdx(:),[3 2 1]);                
                
                T = blkdiag(t,t);
                dTdx = [dtdx zeros(size(dtdx));
                        zeros(size(dtdx)) dtdx];
            
                % Deformations caused by distributed loads
                tempx = 1/(120*E)*NeDLoad*diag([0 0 0 0 0 0])*T;
                tempy = 1/(120*E)*NeDLoad*diag([0 1/(Izz+1e-250) 0 0 1/(Izz+1e-250) 0])*T;
                tempz = 1/(120*E)*NeDLoad*diag([0 0 1/(Iyy+1e-250) 0 0 1/(Iyy+1e-250)])*T;
                %
                Cx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),((irow-1)*6+1):(irow*6)) = t(1,1)*tempx+t(2,1)*tempy+t(3,1)*tempz;
                Cy(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),((irow-1)*6+1):(irow*6)) = t(1,2)*tempx+t(2,2)*tempy+t(3,2)*tempz;
                Cz(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),((irow-1)*6+1):(irow*6)) = t(1,3)*tempx+t(2,3)*tempy+t(3,3)*tempz;
                
                dtempxdx = zeros(nPoints,6);
                for n=1:nVar
                    if any(reshape(dtdx(:,:,n),[],1)~=0)
                        dtempydx = 1/(120*E)*(dNeDLoaddx(:,:,n)*diag([0 1/(Izz+1e-250) 0 0 1/(Izz+1e-250) 0])*T + NeDLoad*diag([0 -dIzzdx(n)/(Izz^2+1e-250) 0 0 -dIzzdx(n)/(Izz^2+1e-250) 0])*T + NeDLoad*diag([0 1/(Izz+1e-250) 0 0 1/(Izz+1e-250) 0])*dTdx(:,:,n));
                        dtempzdx = 1/(120*E)*(dNeDLoaddx(:,:,n)*diag([0 0 1/(Iyy+1e-250) 0 0 1/(Iyy+1e-250)])*T + NeDLoad*diag([0 0 -dIyydx(n)/(Iyy^2+1e-250) 0 0 -dIyydx(n)/(Iyy^2+1e-250)])*T + NeDLoad*diag([0 0 1/(Iyy+1e-250) 0 0 1/(Iyy+1e-250)])*dTdx(:,:,n));
                        %
                        dCxdx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),((irow-1)*6+1):(irow*6),n) = t(1,1)*dtempxdx+t(2,1)*dtempydx+t(3,1)*dtempzdx + dtdx(1,1,n)*tempx+dtdx(2,1,n)*tempy+dtdx(3,1,n)*tempz;    
                        dCydx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),((irow-1)*6+1):(irow*6),n) = t(1,2)*dtempxdx+t(2,2)*dtempydx+t(3,2)*dtempzdx + dtdx(1,2,n)*tempx+dtdx(2,2,n)*tempy+dtdx(3,2,n)*tempz;
                        dCzdx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),((irow-1)*6+1):(irow*6),n) = t(1,3)*dtempxdx+t(2,3)*dtempydx+t(3,3)*dtempzdx + dtdx(1,3,n)*tempx+dtdx(2,3,n)*tempy+dtdx(3,3,n)*tempz;  
                    else
                        dtempydx = 1/(120*E)*(dNeDLoaddx(:,:,n)*diag([0 1/(Izz+1e-250) 0 0 1/(Izz+1e-250) 0])*T + NeDLoad*diag([0 -dIzzdx(n)/(Izz^2+1e-250) 0 0 -dIzzdx(n)/(Izz^2+1e-250) 0])*T);
                        dtempzdx = 1/(120*E)*(dNeDLoaddx(:,:,n)*diag([0 0 1/(Iyy+1e-250) 0 0 1/(Iyy+1e-250)])*T + NeDLoad*diag([0 0 -dIyydx(n)/(Iyy^2+1e-250) 0 0 -dIyydx(n)/(Iyy^2+1e-250)])*T);
                        %
                        dCxdx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),((irow-1)*6+1):(irow*6),n) = t(1,1)*dtempxdx+t(2,1)*dtempydx+t(3,1)*dtempzdx;
                        dCydx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),((irow-1)*6+1):(irow*6),n) = t(1,2)*dtempxdx+t(2,2)*dtempydx+t(3,2)*dtempzdx;
                        dCzdx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),((irow-1)*6+1):(irow*6),n) = t(1,3)*dtempxdx+t(2,3)*dtempydx+t(3,3)*dtempzdx; 
                    end
                end
            end
        end 
    end
    
    % Deformations caused by nodal displacements [calculate Ax Ay Az B]  
    T=blkdiag(t,t,t,t);
    tempx=NeLCS*diag([1 0 0 0 0 0 1 0 0 0 0 0])*T;
    tempy=NeLCS*diag([0 1 0 0 0 L 0 1 0 0 0 L])*T;
    tempz=NeLCS*diag([0 0 1 0 L 0 0 0 1 0 L 0])*T;
    %
    Ax(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF)=(t(1,1)*tempx+t(2,1)*tempy+t(3,1)*tempz)*C;
    Ay(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF)=(t(1,2)*tempx+t(2,2)*tempy+t(3,2)*tempz)*C;
    Az(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF)=(t(1,3)*tempx+t(2,3)*tempy+t(3,3)*tempz)*C;
    
    B(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:3)=[NeLCS(:,1), NeLCS(:,7)]*Node(1:2,:);
    B((nPoints+1)*iElem,1)=NaN;         
    
    dTdx = [dtdx                zeros(size(dtdx)) 	zeros(size(dtdx))   zeros(size(dtdx));
            zeros(size(dtdx))   dtdx                zeros(size(dtdx))   zeros(size(dtdx));
            zeros(size(dtdx))   zeros(size(dtdx))   dtdx                zeros(size(dtdx));
            zeros(size(dtdx))   zeros(size(dtdx))   zeros(size(dtdx))   dtdx];
    for n=1:nVar
        if any(reshape(dtdx(:,:,n),[],1)~=0)
            dtempxdx = NeLCS*diag([1 0 0 0 0 0 1 0 0 0 0 0])*dTdx(:,:,n);
            dtempydx = NeLCS*(diag([0 0 0 0 0 dLdx(n) 0 0 0 0 0 dLdx(n)])*T + diag([0 1 0 0 0 L 0 1 0 0 0 L])*dTdx(:,:,n));
            dtempzdx = NeLCS*(diag([0 0 0 0 dLdx(n) 0 0 0 0 0 dLdx(n) 0])*T + diag([0 0 1 0 L 0 0 0 1 0 L 0])*dTdx(:,:,n));
            %
            dAxdx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF,n) = (t(1,1)*dtempxdx+t(2,1)*dtempydx+t(3,1)*dtempzdx + dtdx(1,1,n)*tempx+dtdx(2,1,n)*tempy+dtdx(3,1,n)*tempz)*C;
            dAydx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF,n) = (t(1,2)*dtempxdx+t(2,2)*dtempydx+t(3,2)*dtempzdx + dtdx(1,2,n)*tempx+dtdx(2,2,n)*tempy+dtdx(3,2,n)*tempz)*C;
            dAzdx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF,n) = (t(1,3)*dtempxdx+t(2,3)*dtempydx+t(3,3)*dtempzdx + dtdx(1,3,n)*tempx+dtdx(2,3,n)*tempy+dtdx(3,3,n)*tempz)*C;
        else
            dtempxdx = NeLCS*diag([0 0 0 0 0 0 0 0 0 0 0 0])*T;
            dtempydx = NeLCS*diag([0 0 0 0 0 dLdx(n) 0 0 0 0 0 dLdx(n)])*T;
            dtempzdx = NeLCS*diag([0 0 0 0 dLdx(n) 0 0 0 0 0 dLdx(n) 0])*T;
            %
            dAxdx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF,n) = (t(1,1)*dtempxdx+t(2,1)*dtempydx+t(3,1)*dtempzdx)*C;
            dAydx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF,n) = (t(1,2)*dtempxdx+t(2,2)*dtempydx+t(3,2)*dtempzdx)*C;
            dAzdx(((nPoints+1)*(iElem-1)+1):((nPoints+1)*iElem-1),1:nDOF,n) = (t(1,3)*dtempxdx+t(2,3)*dtempydx+t(3,3)*dtempzdx)*C;
        end
    end
    
end
