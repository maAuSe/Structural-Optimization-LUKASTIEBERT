function [ForcesLCS,ForcesGCS,dForcesLCSdx,dForcesGCSdx]=elemforces(Nodes,Elements,Types,Sections,Materials,DOF,U,DLoads,TLoads,dNodesdx,dSectionsdx,dUdx,dDLoadsdx)

%ELEMFORCES   Compute the element forces.
%
%   [ForcesLCS,ForcesGCS] = ELEMFORCES(Nodes,Elements,Types,Sections,Materials,DOF,U,DLoads,TLoads)
%   [ForcesLCS,ForcesGCS] = ELEMFORCES(Nodes,Elements,Types,Sections,Materials,DOF,U)
%       computes the element forces in the local (beam convention) and the  
%       global (algebraic convention) coordinate system.
%
%   [ForcesLCS,ForcesGCS,dForcesLCSdx,dForcesGCSdx]
%              = ELEMFORCES(Nodes,Elements,Types,Sections,Materials,DOF,U,DLoads,TLoads,
%                                                              dNodesdx,dSectionsdx,dUdx,dDLoadsdx)
%       additionally computes the derivatives of the element forces with
%       respect to the design variables x.
%
%   Nodes           Node definitions          [NodID x y z]
%   Elements        Element definitions       [EltID TypID SecID MatID n1 n2 ...]
%   Types           Element type definitions  {TypID EltName Option1 ...}
%   Sections        Section definitions       [SecID SecProp1 SecProp2 ...]
%   Materials       Material definitions      [MatID MatProp1 MatProp2 ...]
%   DOF             Degrees of freedom        (nDOF * 1)
%   U               Displacements             (nDOF * nLC)
%   DLoads          Distributed loads         [EltID n1globalX n1globalY n1globalZ ...]
%   TLoads          Temperature loads         [EltID Tyt Tyb Tzt Tzb]
%   dNodesdx        Node definitions derivatives    (SIZE(Nodes) * nVar)
%   dSectionsdx     Section definitions derivatives (SIZE(Sections) * nVar)
%   dUdx            Displacements derivatives       (nDOF * nLC * nVar)
%   dDLoadsdx       Distributed loads derivative	(SIZE(DLoads) * nVar)
%   ForcesLCS       Element forces in LCS (beam convention)  [N Vy Vz T My Mz] (nElem * 12 * nLC)
%   ForcesGCS       Element forces in GCS (algebraic convention)               (nElem * 12 * nLC)
%   dForcesLCSdx    Element forces derivatives in LCS  (nElem * 12 * nLC * nVar)
%   dForcesGCSdx    Element forces derivatives in GCS  (nElem * 12 * nLC * nVar)
%
%   See also FORCES_TRUSS, FORCES_BEAM.

% David Dooms, Wouter Dillen
% October 2008, April 2020

if nargin<8, DLoads = []; end
if nargin<9, TLoads = []; end
if nargin<10, dNodesdx = []; end
if nargin<11, dSectionsdx = []; end
if nargin<12, dUdx = []; end
if nargin<13, dDLoadsdx = []; end

nVar = 0;
if nargout>2 && (~isempty(dNodesdx) || ~isempty(dSectionsdx) || ~isempty(dUdx) || ~isempty(dDLoadsdx))
    nVar = max([size(dNodesdx,3),size(dSectionsdx,3),size(dUdx,3),size(dDLoadsdx,4)]);
end
  
if ~isempty(dNodesdx) && ~isempty(dSectionsdx) && (size(dNodesdx,3)~=size(dSectionsdx,3))
    error('SIZE(dNodesdx,3) and SIZE(dSectionsdx,3) must be identical.');
end

if nVar==0 || isempty(dNodesdx), dNodesdx = zeros([size(Nodes),nVar]); end
if nVar==0 || isempty(dSectionsdx), dSectionsdx = zeros([size(Sections),nVar]); end
if nVar==0 || isempty(dUdx), dUdx = zeros(size(U,1),size(U,2),nVar); end
if nVar==0 || isempty(dDLoadsdx), dDLoadsdx = zeros(size(DLoads,1),size(DLoads,2),size(DLoads,3),nVar); end


nElem=size(Elements,1);
nTimeSteps=size(U,2);

ForcesLCS=zeros(nElem,12,nTimeSteps);
ForcesGCS=zeros(nElem,12,nTimeSteps);
dForcesLCSdx=zeros(nElem,12,nTimeSteps,nVar);
dForcesGCSdx=zeros(nElem,12,nTimeSteps,nVar);


for iElem=1:nElem

    % Type
    TypID=Elements(iElem,2);
    loc=find(cell2mat(Types(:,1))==TypID);
    if isempty(loc)
        error('Element type %i is not defined.',TypID)
    elseif length(loc)>1
        error('Element type %i is multiply defined.',TypID)
    end
    
    Type=Types{loc,2};
    
    if size(Types,2)<3
        Options={};
    else
        Options=Types{loc,3};
    end
    
    % Section
    SecID=Elements(iElem,3);
    loc=find(Sections(:,1)==SecID);
    if isempty(loc)
        error('Section %i is not defined.',SecID)
    elseif length(loc)>1
        error('Section %i is multiply defined.',SecID)
    end
    Section=Sections(loc,2:end);
    dSectiondx=dSectionsdx(loc,2:end,:);
    
    % Material
    MatID=Elements(iElem,4);
    loc=find(Materials(:,1)==MatID);
    if isempty(loc)
        error('Material %i is not defined.',MatID)
    elseif length(loc)>1
        error('Material %i is multiply defined.',MatID)
    end
    Material=Materials(loc,2:end);
    
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
    
    % Displacements
    [UeGCS,LL]=elemdisp(Type,NodeNum,DOF,U);
    dUeGCSdx=zeros(size(UeGCS,1),size(UeGCS,2),nVar);
    for n=1:nVar
        dUeGCSdx(:,:,n)=LL*dUdx(:,:,n);
    end
    
    EltID=Elements(iElem,1);
    
    % DLoads
    if isempty(DLoads)
        DLoads=zeros(1,7,nTimeSteps);
    end
    loc=find(DLoads(:,1)==EltID);
    if isempty(loc)
        DLoad=zeros(6,nTimeSteps);
        dDLoaddx=zeros(6,nTimeSteps,nVar);
    else
        DLoad=permute(DLoads(loc,2:end,:),[2 3 1]);
        dDLoaddx=permute(dDLoadsdx(loc,2:end,:,:),[2 3 4 1]);
    end
    
    % TLoads
    if isempty(TLoads)
        TLoads=zeros(1,4,nTimeSteps);
    end
    
    % Define neutral line to compute average temperature (in case of
    % asymmetric profiles)
    if length(Section)>8
        yt=Section(7); yb=Section(8);
        zt=Section(9); zb=Section(10);
    elseif length(Section)>6
        yt=Section(7); yb=Section(8);
        zt=1; zb=1;
    else
        yt=1; yb=1;
        zt=1; zb=1;
    end
    if strcmp(Type,'truss')
        yt=1; yb=1;
        zt=1; zb=1;
    end
    
    loc=find(TLoads(:,1)==EltID);
    if isempty(loc)
        TLoad=zeros(3,nTimeSteps);
    elseif length(loc)>1
        error('Element %i has multiple temperature loads.',EltID)
    else
        nTemp=size(TLoads(loc,2:end,:),2);
        TLoad=zeros(1,3,nTimeSteps);
        if nTemp==1
            TLoad(1,1,:)=TLoads(loc,2,:);
        elseif nTemp==2
            % Average temperature 
            Tmy=yb/(yt+yb)*TLoads(loc,2,:)+yt/(yt+yb)*TLoads(loc,3,:);
            TLoad(1,1,:)=Tmy;
            % Temperature gradient Tyt-tyb (local y direction only)
            TLoad(1,2,:)=TLoads(loc,2,:)-TLoads(loc,3,:);
        elseif nTemp==4
            % Check for 2 zero temperatures in either the local y or z direction
            % since this influences the average temperature: i0 = 0, 1 or 2
            i0=any(TLoads(loc,2:3,:),2)+any(TLoads(loc,4:5,:),2);
            % Weight function for average temperatures in local y or z direction: 
            % projects 0, 1 or 2 on 0, 1 or 1/2 respectively
            w=(-3/4*(i0).^2+7/4*i0);    
            % Average temperature
            Tmy=yb/(yt+yb)*TLoads(loc,2,:)+yt/(yt+yb)*TLoads(loc,3,:);
            Tmz=zb/(zt+zb)*TLoads(loc,4,:)+zt/(zt+zb)*TLoads(loc,5,:);
            TLoad(1,1,:)=w.*(Tmy+Tmz);
            % Temperature gradient (local y and z direction)
            TLoad(1,2,:)=TLoads(loc,2,:)-TLoads(loc,3,:); % y dir
            TLoad(1,3,:)=TLoads(loc,4,:)-TLoads(loc,5,:); % z dir
        else
            error('Check temperature definitions.')
        end
        TLoad=permute(TLoad,[2 3 1]);
    end
           
    % Forces
    [temp1,temp2,temp3,temp4] = eval(['forces_' Type '(Node,Section,Material,UeGCS,DLoad,TLoad,Options,dNodedx,dSectiondx,dUeGCSdx,dDLoaddx)']);
    ForcesLCS(iElem,:,:) = permute(temp1,[3 1 2]);
    ForcesGCS(iElem,:,:) = permute(temp2,[3 1 2]);
    dForcesLCSdx(iElem,:,:,:) = permute(temp3,[4 1 2 3]);
    dForcesGCSdx(iElem,:,:,:) = permute(temp4,[4 1 2 3]);
            
end

% change to beam convention
ForcesLCS(:,[1:5,12],:) = -ForcesLCS(:,[1:5,12],:);
dForcesLCSdx(:,[1:5,12],:,:) = -dForcesLCSdx(:,[1:5,12],:,:);
    
