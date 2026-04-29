function [UxdiagrGCS,UydiagrGCS,UzdiagrGCS,dUxdiagrGCSdx,dUydiagrGCSdx,dUzdiagrGCSdx] = udiagrgcs(Nodes,Elements,Types,DOF,U,DLoads,Sections,Materials,Points,dNodesdx,dUdx,dDLoadsdx,dSectionsdx) 

%UDIAGRGCS   Return displacement diagrams in GCS
%
%   [UxdiagrGCS,UydiagrGCS,UzdiagrGCS] = UDIAGRGCS(Nodes,Elements,Types,DOF,U,DLoads,
%                                                               Sections,Materials,Points)
%   [UxdiagrGCS,UydiagrGCS,UzdiagrGCS] = UDIAGRGCS(Nodes,Elements,Types,DOF,U,DLoads,
%                                                               Sections,Materials)
%   [UxdiagrGCS,UydiagrGCS,UzdiagrGCS] = UDIAGRGCS(Nodes,Elements,Types,DOF,U,[],
%                                                               Sections,Materials)
%   [UxdiagrGCS,UydiagrGCS,UzdiagrGCS] = UDIAGRGCS(Nodes,Elements,Types,DOF,U)
%       computes the displacements of the interpolation points after
%       deformation in the global (algebraic) coordinate system. If DLoads,
%       Sections and Materials are supplied, the displacements that occur
%       due to distributed loads if all nodes are fixed, are superimposed.
%
%   [UxdiagrGCS,UydiagrGCS,UzdiagrGCS,dUxdiagrGCSdx,dUydiagrGCSdx,dUzdiagrGCSdx] 
%           = UDIAGRGCS(Nodes,Elements,Types,DOF,U,DLoads,Sections,Materials,Points,
%                                                     dNodesdx,dUdx,dDLoadsdx,dSectionsdx)
%       additionally computes the derivatives of the displacement values
%       with respect to the design variables x.
%
%   Nodes      Node definitions          [NodID x y z]
%   Elements   Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   DOF        Degrees of freedom  (nDOF * 1)
%   U          Displacements  (nDOF * 1)
%   DLoads     Distributed loads         [EltID n1globalX n1globalY n1globalZ ...]
%                   (use empty array [] when shear deformation (in beam element) 
%                   is considered but no DLoads are present)
%   Sections   Section definitions       [SecID SecProp1 SecProp2 ...]
%   Materials  Material definitions      [MatID MatProp1 MatProp2 ... ]
%   Points     Points in the local coordinate system            (1 * nPoints)
%   dNodesdx        Node definitions derivatives                (SIZE(Nodes) * nVar)
%   dUdx            Displacements derivatives                   (SIZE(U) * nVar)
%   dDLoadsdx       Distributed loads derivatives               (SIZE(DLoads) * nVar)
%   UxdiagrGCS x-direction displacement values at the points    (nElem * nPoints * nLC)
%   UydiagrGCS y-direction displacement values at the points    (nElem * nPoints * nLC)
%   UzdiagrGCS z-direction displacement values at the points    (nElem * nPoints * nLC)
%   dUxdiagrGCSdx   x-direction displacement values derivatives	(nElem * nPoints * nLC * nVar)
%   dUydiagrGCSdx   y-direction displacement values derivatives	(nElem * nPoints * nLC * nVar)
%   dUzdiagrGCSdx   z-direction displacement values derivatives	(nElem * nPoints * nLC * nVar)
%
%   See also PLOTDISP, DISP_TRUSS, DISP_BEAM.

% Wouter Dillen
% July 2017

if nargin<6, DLoads = []; end
if nargin<7, Sections = []; end
if nargin<8, Materials = []; end
if nargin<9, Points = []; end
if nargin<10, dNodesdx = []; end
if nargin<11, dUdx = []; end
if nargin<12, dDLoadsdx = []; end
if nargin<13, dSectionsdx = []; end

nVar = 0;
if nargout>3 && (~isempty(dNodesdx) || ~isempty(dUdx) || ~isempty(dDLoadsdx) || ~isempty(dSectionsdx))
    nVar = max([size(dNodesdx,3),size(dUdx,3),size(dDLoadsdx,4),size(dSectionsdx,3)]);
end

if isempty(Points)
    nPoints=21;
    Points=linspace(0,1,nPoints);
else
    nPoints=length(Points);
end
if nVar==0 || isempty(dNodesdx), dNodesdx = zeros([size(Nodes),nVar]); end
if nVar==0 || isempty(dUdx), dUdx = zeros(size(U,1),size(U,2),nVar); end
if nVar==0 || isempty(dDLoadsdx), dDLoadsdx = zeros(size(DLoads,1),size(DLoads,2),size(DLoads,3),nVar); end
if nVar==0 || isempty(dSectionsdx), dSectionsdx = zeros([size(Sections),nVar]); end


nType = size(Types,1);
nElem = size(Elements,1);
nLC = size(U,2);

UxdiagrGCS = zeros(nElem,nPoints,nLC);
UydiagrGCS = zeros(nElem,nPoints,nLC);
UzdiagrGCS = zeros(nElem,nPoints,nLC);
dUxdiagrGCSdx = zeros(nElem,nPoints,nLC,nVar);
dUydiagrGCSdx = zeros(nElem,nPoints,nLC,nVar);
dUzdiagrGCSdx = zeros(nElem,nPoints,nLC,nVar);

Ax=[];
Ay=[];
Az=[];
B=[];
dAxdx=[];
dAydx=[];
dAzdx=[];

if ~isempty(DLoads)
    Cx=[];
    Cy=[];
    Cz=[];
    dCxdx=[];
    dCydx=[];
    dCzdx=[];
end

ElemOrder = [];

% Compute the interpolation point coordinates of the deformed elements
for iType=1:nType
    loc=find(Elements(:,2)==cell2mat(Types(iType,1)));
    if isempty(loc)
        
    else
        ElemType=Elements(loc,:);
        Type=Types{iType,2};
        ElemOrder=[ElemOrder; Elements(loc,1)];
        
        if ~isempty(DLoads)
            [Axadd,Ayadd,Azadd,Badd,Cxadd,Cyadd,Czadd,dAxadddx,dAyadddx,dAzadddx,dCxadddx,dCyadddx,dCzadddx] = ...
                eval(['disp_' Type '(Nodes,ElemType,DOF,DLoads,Sections,Materials,Points,dNodesdx,dDLoadsdx,dSectionsdx)']);
            Cx=[Cx; Cxadd];
            Cy=[Cy; Cyadd];
            Cz=[Cz; Czadd];
            dCxdx=cat(1,dCxdx,dCxadddx);
            dCydx=cat(1,dCydx,dCyadddx);
            dCzdx=cat(1,dCzdx,dCzadddx);
        else
            [Axadd,Ayadd,Azadd,Badd,~,~,~,dAxadddx,dAyadddx,dAzadddx] = ...
                eval(['disp_' Type '(Nodes,ElemType,DOF,[],Sections,Materials,Points,dNodesdx,[],dSectionsdx)']);
        end 
        Ax=[Ax; Axadd];
        Ay=[Ay; Ayadd];
        Az=[Az; Azadd];
        B=[B; Badd];
        dAxdx=cat(1,dAxdx,dAxadddx);
        dAydx=cat(1,dAydx,dAyadddx);
        dAzdx=cat(1,dAzdx,dAzadddx);
    end
end

if ~isempty(DLoads)
    DLoads=permute(DLoads(:,2:7,:),[2 1 3]);
    DLoad=reshape(DLoads,numel(DLoads(:,:,1)),size(DLoads,3));
    X=Ax*U+Cx*DLoad;
    Y=Ay*U+Cy*DLoad;
    Z=Az*U+Cz*DLoad;
    
    dDLoadsdx=permute(dDLoadsdx(:,2:7,:,:),[2 1 3 4]);
    dDLoaddx=reshape(dDLoadsdx,size(dDLoadsdx,1)*size(dDLoadsdx,2),size(dDLoadsdx,3),size(dDLoadsdx,4));
    dXdx=zeros([size(X),nVar]); dYdx=zeros([size(Y),nVar]); dZdx=zeros([size(Z),nVar]);
    for n=1:nVar
        dXdx(:,:,n)=dAxdx(:,:,n)*U+Ax*dUdx(:,:,n)+dCxdx(:,:,n)*DLoad+Cx*dDLoaddx(:,:,n);
        dYdx(:,:,n)=dAydx(:,:,n)*U+Ay*dUdx(:,:,n)+dCydx(:,:,n)*DLoad+Cy*dDLoaddx(:,:,n);
        dZdx(:,:,n)=dAzdx(:,:,n)*U+Az*dUdx(:,:,n)+dCzdx(:,:,n)*DLoad+Cz*dDLoaddx(:,:,n);
    end
else
    X=Ax*U;
    Y=Ay*U;
    Z=Az*U;
    
    dXdx=zeros([size(X),nVar]); dYdx=zeros([size(Y),nVar]); dZdx=zeros([size(Z),nVar]);
    for n=1:nVar
        dXdx(:,:,n)=dAxdx(:,:,n)*U+Ax*dUdx(:,:,n);
        dYdx(:,:,n)=dAydx(:,:,n)*U+Ay*dUdx(:,:,n);
        dZdx(:,:,n)=dAzdx(:,:,n)*U+Az*dUdx(:,:,n);
    end
end

% Reformulate coordinate information to comply with the output matrix structure
[~,order] = ismember(ElemOrder,Elements(:,1));
for iElem=1:length(ElemOrder)
    indices = (iElem-1)*(nPoints+1)+1 : (iElem-1)*(nPoints+1)+nPoints;
    UxdiagrGCS(order(iElem),:,:) = permute(X(indices,:),[3 1 2]);
    UydiagrGCS(order(iElem),:,:) = permute(Y(indices,:),[3 1 2]);
    UzdiagrGCS(order(iElem),:,:) = permute(Z(indices,:),[3 1 2]);
    dUxdiagrGCSdx(order(iElem),:,:,:) = permute(dXdx(indices,:,:),[4 1 2 3]);
    dUydiagrGCSdx(order(iElem),:,:,:) = permute(dYdx(indices,:,:),[4 1 2 3]);
    dUzdiagrGCSdx(order(iElem),:,:,:) = permute(dZdx(indices,:,:),[4 1 2 3]);
end

end