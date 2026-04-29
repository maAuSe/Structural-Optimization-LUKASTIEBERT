function [SeGCS,SeLCS,vLCS] = se_shell4(Node,Section,Material,UeGCS,Options,gcs)

%SE_SHELL4   Compute the element stresses for a shell4 element.
%
%   [SeGCS,SeLCS,vLCS] = se_shell4(Node,Section,Material,UeGCS,Options,GCS)
%   [SeGCS,SeLCS]      = se_shell4(Node,Section,Material,UeGCS,Options,GCS)
%    SeGCS             = se_shell4(Node,Section,Material,UeGCS,Options,GCS)
%   computes the element stresses in the global and the  
%   local coordinate system for the shell4 element.
%
%   Node       Node definitions           [x y z] (4 * 3)
%              Nodes should have the following order:
%   Section    Section definition         [h]
%   Material   Material definition        [E nu rho]
%   UeGCS      Displacements (24 * nTimeSteps)
%   Options    Element options            {Option1 Option2 ...}
%   GCS        Global coordinate system in which stresses are returned
%              'cart'|'cyl'|'sph'
%   SeGCS      Element stresses in GCS in corner nodes IJKL and 
%              at top/mid/bot of shell (72 * nTimeSteps)   
%              72 = 6 stress comp. * 4 nodes * 3 locations
%                                        [sxx syy szz sxy syz sxz]
%   SeLCS      Element stresses in LCS in corner nodes IJKL and 
%              at top/mid/bot of shell (72 * nTimeSteps)   
%                                        [sxx syy szz sxy syz sxz]
%   vLCS       Unit vectors of LCS (1 * 9)
%
%   See also ELEMSTRESS, SE_SHELL8.

% Miche Jansen
% 2010

% Material
E=Material(1,1);
nu=Material(1,2);

% Section
h=Section(1,1);

% Check nodes
Node=Node(1:4,1:3);
if ~ all(isfinite(Node))
    error('Not all the nodes exist.')
end

% Transformation matrix
[t,Node_lc,W]=trans_shell4(Node);
T=blkdiag(t,t,t,t,t,t,t,t);

UeLCS=W.'*T*UeGCS;

SeLCS = selcs_shell4(Node_lc,h,E,nu,UeLCS,Options);

% transformatie matrix uit p194 Zienkiewicz deel1 (Tsigma^(-1)=Tespilon^(T))
theta = vtrans_solid(t.'); % transformation LCS -> GCS
   
SeGCS = zeros(size(SeLCS));
for ind = 1:12; 
    SeGCS((1:6)+6*(ind-1),:)= theta*SeLCS((1:6)+6*(ind-1),:);
end

% LCS
if nargout > 2
vLCS = t.';
vLCS = vLCS(:).';
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