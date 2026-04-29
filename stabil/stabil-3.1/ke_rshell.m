function [Ke,Me,dKedx] = ke_rshell(Node,Section,Material,Options,dNodedx,dSectiondx)

% Mattias Schevenels
% April 2008

if nargin<5, dNodedx = []; end
if nargin<6, dSectiondx = []; end
if (nnz(dNodedx)+nnz(dSectiondx)) > 0
    error('Sensitivities have not been implemented yet.')
end
nVar = 0;
if nargout>2 && (~isempty(dNodedx) || ~isempty(dSectiondx))
    nVar = max(size(dNodedx,3),size(dSectiondx,3));
end

if nargout>1, Me = []; end
if nargout>2, dKedx = cell(nVar,1); end


% Check nodes
if ~ all(isfinite(Node(1:4,1:3)))
    error('Not all the nodes exist.')
end

% CHECK ELEMENT SHAPE
Lc=max(Node(:))-min(Node(:));
if Lc==0
  error('RSHELL elements must be rectangular.');
end
if abs(dot(Node(2,:)-Node(1,:),Node(4,:)-Node(1,:)))/Lc>1e-10
  error('RSHELL elements must be rectangular.');
end
if norm(Node(2,:)-Node(1,:)-Node(3,:)+Node(4,:))/Lc>1e-10
  error('RSHELL elements must be rectangular.');
end

% ELEMENT SIZE
Lx=norm(Node(2,:)-Node(1,:));
Ly=norm(Node(4,:)-Node(1,:));

% MATERIAL PROPERTIES
E=Material(1);
nu=Material(2);

% SECTION PROPERTIES
t=Section(1);
d=Section(2);

% TRANSFORMATION MATRIX
Tn=trans_beam(Node);
T=blkdiag(Tn,Tn,Tn,Tn,Tn,Tn,Tn,Tn);

if nargout==2
  % COMPUTE STIFFNESS AND MASS MATRICES
  rho=Material(3);
  [KeLCS,MeLCS]=kelcs_rshell(Lx,Ly,t,d,E,nu,rho);
  Me=T.'*MeLCS*T;
  Ke=T.'*KeLCS*T;
else
  % COMPUTE STIFFNESS MATRIX
  KeLCS=kelcs_rshell(Lx,Ly,t,d,E,nu);
  Ke=T.'*KeLCS*T;
end

