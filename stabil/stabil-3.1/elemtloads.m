function F=elemtloads(TLoads,Nodes,Elements,Types,Sections,Materials,DOF)

%ELEMTLOADS   Equivalent nodal forces for temperature loading.
%
%   F=elemtloads(TLoads,Nodes,Elements,Types,Sections,Materials,DOF)
%   computes the equivalent nodal forces of a temperature gradient
%   (in the global coordinate system).
%
%   TLoads     Temperature gradient         [EltID Tyt Tyb Tzt Tzb]
%              Tkt and Tkb correspond to the temperatures at the top and
%              the bottom of the profile when the k-axis (LCS) points up (k=y,z)
%   Nodes      Node definitions             [NodID x y z]
%   Elements   Element definitions          [EltID TypID SecID MatID n1 n2 ...]
%   Types      Type definitions             {TypID EltName Option1 ... }
%   Sections   Section definitions          [SecID A ky kz Ixx Iyy Izz yt yb zt zb]
%   Materials  Material definitions         [MatID E nu rho alpha]
%   DOF        Degrees of freedom           (nDOF * 1)
%   F          Load vector                  (nDOF * 1)
%
%   See also TLOADS_TRUSS, TLOADS_BEAM.

% Pieter Reumers, Jef Wambacq
% 2018

nTLoads=size(TLoads,1);

nterm=0;

for iTLoad=1:nTLoads
    EltID=TLoads(iTLoad,1);
    iElem=find(Elements(:,1)==EltID);
    
    Type=Types{Elements(iElem,2),2};
    Section=Sections(Elements(iElem,3),2:end);
    Material=Materials(Elements(iElem,4),2:end);
    NodeNum=Elements(iElem,5:end);
    Node=zeros(length(NodeNum),3);
    for iNode=1:length(NodeNum)
        loc=find(Nodes(:,1)==NodeNum(1,iNode));
        if isempty(loc)
            Node(iNode,:)=[NaN NaN NaN];
        elseif length(loc)>1
            error('Node %i is multiply defined.',NodeNum(1,iNode))
        else
            Node(iNode,:)=Nodes(loc,2:end);
        end
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
    
    nTemp=size(TLoads(iTLoad,2:end,:),2);   % Number of input temperatures
    nLC=size(TLoads(iTLoad,2:end,:),3);     % Number of load steps/cases
    TLoad=zeros(1,3,nLC);
    if nTemp==1
        TLoad(1,1,:)=TLoads(iTLoad,2,:);
    elseif nTemp==2
        % Average temperature 
        Tmy=yb/(yt+yb)*TLoads(iTLoad,2,:)+yt/(yt+yb)*TLoads(iTLoad,3,:);
        TLoad(1,1,:)=Tmy;
        % Temperature gradient (local y direction only)
        TLoad(1,2,:)=TLoads(iTLoad,2,:)-TLoads(iTLoad,3,:);
    elseif nTemp==4
        % Check for 2 zero temperatures in either the local y or z direction
        % since this influences the average temperature: i0 = 0, 1 or 2
        i0=any(TLoads(iTLoad,2:3,:),2)+any(TLoads(iTLoad,4:5,:),2);
        % Weight function for average temperatures in local y or z direction: 
        % projects 0, 1 or 2 on 0, 1 or 1/2 respectively
        w=(-3/4*(i0).^2+7/4*i0);    
        % Average temperature
        Tmy=yb/(yt+yb)*TLoads(iTLoad,2,:)+yt/(yt+yb)*TLoads(iTLoad,3,:);
        Tmz=zb/(zt+zb)*TLoads(iTLoad,4,:)+zt/(zt+zb)*TLoads(iTLoad,5,:);
        TLoad(1,1,:)=w.*(Tmy+Tmz);
        % Temperature gradient (local y and z direction)
        TLoad(1,2,:)=TLoads(iTLoad,2,:)-TLoads(iTLoad,3,:); % y dir
        TLoad(1,3,:)=TLoads(iTLoad,4,:)-TLoads(iTLoad,5,:); % z dir
    else
        error('Check temperature definitions.')
    end
    TLoad=permute(TLoad,[2 3 1]);
    
    Fe=eval(['tloads_' Type '(TLoad,Node,Section,Material)']);
    
    PLoad(nterm+1:nterm+size(Fe,1),:)=Fe;
    seldof(nterm+1:nterm+size(Fe,1),1)=eval(['dof_' Type '(NodeNum)']);
    
    nterm=nterm+size(Fe,1);
    
end

F=nodalvalues(DOF,seldof,PLoad);
