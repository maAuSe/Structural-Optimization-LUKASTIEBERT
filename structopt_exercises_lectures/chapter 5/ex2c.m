%EX2A
%
%   Example 5.2: Two-dimensional truss topology design in 88 line style.
%   Solution with GCMMA; sensitivities are only computed when needed.

% Mattias Schevenels
% February 2021

% NODES
[xNode,yNode] = ndgrid(0:8,0:4);
xNode = xNode(:);
yNode = yNode(:);
nNode = length(xNode);
nodes = [xNode yNode];

% ELEMENTS
[node1,node2] = ndgrid(1:nNode,1:nNode);
node1 = triu(node1,1);
node2 = triu(node2,1);
node1 = node1(find(node1~=0));
node2 = node2(find(node2~=0));
elements = [node1 node2];
L = sqrt((xNode(node2)-xNode(node1)).^2+(yNode(node2)-yNode(node1)).^2);

% REMOVE ELEMENTS LONGER THAN LMAX
Lmax = 1.5;
iKeep = L<=Lmax;
elements = elements(iKeep,:);
node1 = node1(iKeep);
node2 = node2(iKeep);
L = L(iKeep);
nElem = size(elements,1);

% REMOVE COINCIDING ELEMENTS
thetaMax = 1/180*pi;
[~,iSort] = sort(L);
elements = elements(iSort,:);
node1 = node1(iSort);
node2 = node2(iSort);
L = L(iSort);
iKeep = true(nElem,1);
for iElem = 1:nElem
  for jElem = 1:iElem-1
    if iKeep(jElem)
      if node1(iElem)==node1(jElem)
        v1 = [xNode(node2(iElem))-xNode(node1(iElem)),yNode(node2(iElem))-yNode(node1(iElem))];
        v1 = v1/norm(v1);
        v2 = [xNode(node2(jElem))-xNode(node1(jElem)),yNode(node2(jElem))-yNode(node1(jElem))];
        v2 = v2/norm(v2);
        if abs(acos(dot(v1,v2)))<thetaMax, iKeep(iElem) = false; end
      end
      if node1(iElem)==node2(jElem)
        v1 = [xNode(node2(iElem))-xNode(node1(iElem)),yNode(node2(iElem))-yNode(node1(iElem))];
        v1 = v1/norm(v1);
        v2 = [xNode(node1(jElem))-xNode(node2(jElem)),yNode(node1(jElem))-yNode(node2(jElem))];
        v2 = v2/norm(v2);
        if abs(acos(dot(v1,v2)))<thetaMax, iKeep(iElem) = false; end
      end
      if node2(iElem)==node1(jElem)
        v1 = [xNode(node1(iElem))-xNode(node2(iElem)),yNode(node1(iElem))-yNode(node2(iElem))];
        v1 = v1/norm(v1);
        v2 = [xNode(node2(jElem))-xNode(node1(jElem)),yNode(node2(jElem))-yNode(node1(jElem))];
        v2 = v2/norm(v2);
        if abs(acos(dot(v1,v2)))<thetaMax, iKeep(iElem) = false; end
      end
      if node2(iElem)==node2(jElem)
        v1 = [xNode(node1(iElem))-xNode(node2(iElem)),yNode(node1(iElem))-yNode(node2(iElem))];
        v1 = v1/norm(v1);
        v2 = [xNode(node1(jElem))-xNode(node2(jElem)),yNode(node1(jElem))-yNode(node2(jElem))];
        v2 = v2/norm(v2);
        if abs(acos(dot(v1,v2)))<thetaMax, iKeep(iElem) = false; end
      end
    end
  end
end
elements = elements(iKeep,:);
node1 = node1(iKeep);
node2 = node2(iKeep);
L = L(iKeep);
nElem = size(elements,1);

% YOUNG'S MODULUS
E = 1;

% ELEMENT STIFFNESS MATRICES IN SPARSE FORM
[iK,jK,sK] = deal(zeros(nElem,16));
edof = zeros(nElem,4);
for iElem = 1:nElem
  theta = atan2(yNode(node2(iElem))-yNode(node1(iElem)),xNode(node2(iElem))-xNode(node1(iElem)));
  c = cos(theta);
  s = sin(theta);
  k0 = E/L(iElem)*[c^2  s*c -c^2 -s*c;
                   s*c  s^2 -s*c -s^2;
                  -c^2 -s*c  c^2  s*c;
                  -s*c -s^2  s*c  s^2];
  sK(iElem,:) = k0(:)';
  edof(iElem,:) = [2*node1(iElem)-[1 0] 2*node2(iElem)-[1 0]];
  iK(iElem,:) = reshape(kron(edof(iElem,:),ones(4,1))',[],1);
  jK(iElem,:) = reshape(kron(edof(iElem,:),ones(1,4))',[],1);
end
k1  = sK(:,1);
k2  = sK(:,2);
k3  = sK(:,3);
k4  = sK(:,4);
k6  = sK(:,6);
k7  = sK(:,7);
k8  = sK(:,8);
k11 = sK(:,11);
k12 = sK(:,12);
k16 = sK(:,16);

% BOUNDARY CONDITIONS
ndof = 2*nNode;
xdof = reshape(kron(xNode',[1;1]),[],1);
ydof = reshape(kron(yNode',[1;1]),[],1);
ddof = repmat(['x';'y'],nNode,1);
fixeddofs = find(xdof==0);
alldofs = [1:ndof]';
freedofs = setdiff(alldofs,fixeddofs);
U = zeros(ndof,1);

% LOADS
F = zeros(ndof,1);
F((xdof==8)&(ydof==0)&(ddof=='y')) = -1;

% FORMULATION OF THE OPTIMIZATION PROBLEM
x = repmat(1,nElem,1);            % Initial design
Vmax = L'*x;                      % Maximum allowed volume
xmin = 0.01;                      % Minimum section area
xmax = 20;                        % Maximum section area

% OPTIMIZATION LOOP
iter = 1;                         % Iteration counter
change = inf;                     % (Fictitious) initial design change
tol = xmin/10;                    % Convergence threshold
mmaparams = [];                   % Internal MMA parameters
while change>tol

  % FINITE ELEMENT ANALYSIS
  K = sparse(iK,jK,sK.*repmat(x,1,16));
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);

  % OBJECTIVE FUNCTION AND CONSTRAINT
  C = F'*U;
  V = L'*x;
  f0 = C;
  f = V/Vmax-1;

  % SENSITIVITIES
  if iter==1 || (f0<=subp.f0e && all(f<=subp.fe))
    u1 = U(edof(:,1));
    u2 = U(edof(:,2));
    u3 = U(edof(:,3));
    u4 = U(edof(:,4));
    dCdx = -(k1.*u1.^2+2*k2.*u1.*u2+2.*k3.*u1.*u3+2*k4.*u1.*u4+k6.*u2.^2+2*k7.*u2.*u3+2*k8.*u2.*u4+k11.*u3.^2+2*k12.*u3.*u4+k16.*u4.^2);
    dVdx = L;
    df0dx = dCdx;
    dfdx = dVdx/Vmax;
  end

  % GCMMA UPDATE
  [xnew,y,z,lmult,mmaparams,subp,change,history] = gcmma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams);

  % PLOT CURRENT DESIGN
  if ~isinf(change)
    h = plot([xNode(node1) xNode(node2)]',[yNode(node1) yNode(node2)]','k');
    [~,iSort] = sort(x,1,'descend');
    set(gca,'Children',h(iSort));
    x0 = max(x);
    arrayfun(@(h,x)set(h,'LineWidth',x/x0*8),h,x);
    arrayfun(@(h,x)set(h,'Color',max(0,min(1,(x-x0/16)/(-x0/16)))*[1 1 1]),h,x);
    set(gcf,'Color',[1 1 1]);
    axis('equal');
    axis('off');
    drawnow;
  end

  % CHECK CONVERGENCE CRITERION AND MOVE ON TO NEXT ITERATION
  if change>tol
    iter = iter+1;
    x = xnew;
  end

end
