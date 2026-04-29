function [Nodes,dNodesdx,NodesT,NodesL,NodesB,NodesR] = coonsmesh(NodesT,NodesL,NodesB,NodesR,dNodesTdx,dNodesLdx,dNodesBdx,dNodesRdx)

%COONSMESH   Generate structured Coons surface mesh and sensitivities.
%
%   [Nodes,dNodesdx] = COONSMESH(NodesT,NodesL,NodesB,NodesR,dNodesTdx,dNodesLdx,dNodesBdx,dNodesRdx)
%   generates a structured mesh by transfinite interpolation (Coons surface).  
%   The boundary is provided as four ordered boundary node lists (top/left/bottom/right).  
%   If sensitivities of the boundary nodes with respect to design variables are 
%   provided, the function also returns sensitivities of all interior nodes.
%
%   [Nodes,dNodesdx,NodesT,NodesL,NodesB,NodesR] = COONSMESH(...) additionally 
%   returns the boundary nodes, with their newly assigned IDs in the first column.
%
%   INPUT ARGUMENTS
%
%   NodesT        Nodes on the top boundary, sorted from left to right (Nu * 4).
%   NodesL        Nodes on the left boundary, sorted from bottom to top (Nv * 4).
%   NodesB        Nodes on the bottom boundary, sorted from left to right (Nu * 4).
%   NodesR        Nodes on the right boundary, sorted from bottom to top (Nv * 4).
%                 Corner nodes must coincide consistently between boundaries. 
%                 The first column (node ID) is ignored; only columns 2:4 are used.
%                 
%   dNodesTdx     NodesT derivatives (Nu * 4 * nVar) or [].
%   dNodesLdx     NodesL derivatives (Nv * 4 * nVar) or [].
%   dNodesBdx     NodesB derivatives (Nu * 4 * nVar) or [].
%   dNodesRdx     NodesR derivatives (Nv * 4 * nVar) or [].
%                 If omitted or empty, dNodesdx is returned empty.
%
%   OUTPUT ARGUMENTS
%
%   Nodes         Generated nodes (Nu*Nv * 4).
%                 Nodes are numbered first in the u-direction (left to right), then
%                 in the v-direction (bottom to top).
%   dNodesdx      Nodes derivatives (Nu*Nv * 4 * nVar).
%   NodesT/L/B/R  Same as input arguments but with new node IDs in the first column.

% Mattias Schevenels
% January 2026

% INPUT CHECKS AND DEFAULT VALUES
if any([size(NodesT,2),size(NodesL,2),size(NodesB,2),size(NodesR,2)]~=4)
  error('NodesT/L/B/R must be (N * 4) matrices [id x y z].');
end
Nu = size(NodesT,1);  Nv = size(NodesL,1);
if size(NodesB,1)~=Nu, error('NodesB must have the same number of rows as NodesT.'); end
if size(NodesR,1)~=Nv, error('NodesR must have the same number of rows as NodesL.'); end

t = NodesT(:,2:4);  l = NodesL(:,2:4);  b = NodesB(:,2:4);  r = NodesR(:,2:4);

allxyz = [t; l; b; r];
scale  = max(max(allxyz)-min(allxyz));
tol    = 1e-9*max(1,scale);

if norm(b(1,:)-l(1,:))     > tol, error('Bottom-left corner mismatch: NodesB(1,2:4) vs NodesL(1,2:4).'); end
if norm(b(end,:)-r(1,:))   > tol, error('Bottom-right corner mismatch: NodesB(end,2:4) vs NodesR(1,2:4).'); end
if norm(t(1,:)-l(end,:))   > tol, error('Top-left corner mismatch: NodesT(1,2:4) vs NodesL(end,2:4).'); end
if norm(t(end,:)-r(end,:)) > tol, error('Top-right corner mismatch: NodesT(end,2:4) vs NodesR(end,2:4).'); end

% OPTIONAL SENSITIVITIES: allow [] for boundaries independent of design variables
if nargin < 8, dNodesRdx = []; end
if nargin < 7, dNodesBdx = []; end
if nargin < 6, dNodesLdx = []; end
if nargin < 5, dNodesTdx = []; end

haveSens = ~isempty(dNodesTdx) || ~isempty(dNodesLdx) || ~isempty(dNodesBdx) || ~isempty(dNodesRdx);

if haveSens
  % Determine nVar and Expand empty sensitivities to zeros with correct dimensions
  nVar = max([size(dNodesTdx,3), size(dNodesLdx,3), size(dNodesBdx,3), size(dNodesRdx,3)]);
  if isempty(dNodesTdx), dNodesTdx = zeros(Nu,4,nVar); end
  if isempty(dNodesBdx), dNodesBdx = zeros(Nu,4,nVar); end
  if isempty(dNodesLdx), dNodesLdx = zeros(Nv,4,nVar); end
  if isempty(dNodesRdx), dNodesRdx = zeros(Nv,4,nVar); end

  % Dimension checks
  if any([ndims(dNodesTdx),ndims(dNodesLdx),ndims(dNodesBdx),ndims(dNodesRdx)]~=3)
    error('dNodesXdx must be 3D arrays (N * 4 * nVar).');
  end
  if any([size(dNodesTdx,1),size(dNodesTdx,2),size(dNodesTdx,3)]~=[Nu,4,nVar]), error('dNodesTdx must be (Nu * 4 * nVar).'); end
  if any([size(dNodesBdx,1),size(dNodesBdx,2),size(dNodesBdx,3)]~=[Nu,4,nVar]), error('dNodesBdx must be (Nu * 4 * nVar).'); end
  if any([size(dNodesLdx,1),size(dNodesLdx,2),size(dNodesLdx,3)]~=[Nv,4,nVar]), error('dNodesLdx must be (Nv * 4 * nVar).'); end
  if any([size(dNodesRdx,1),size(dNodesRdx,2),size(dNodesRdx,3)]~=[Nv,4,nVar]), error('dNodesRdx must be (Nv * 4 * nVar).'); end

  % Corner sensitivity consistency
  sBL = reshape(dNodesBdx(1,  2:4,:),3,nVar) - reshape(dNodesLdx(1,  2:4,:),3,nVar);
  sBR = reshape(dNodesBdx(end,2:4,:),3,nVar) - reshape(dNodesRdx(1,  2:4,:),3,nVar);
  sTL = reshape(dNodesTdx(1,  2:4,:),3,nVar) - reshape(dNodesLdx(end,2:4,:),3,nVar);
  sTR = reshape(dNodesTdx(end,2:4,:),3,nVar) - reshape(dNodesRdx(end,2:4,:),3,nVar);
  if max(abs(sBL(:))) > tol, error('Bottom-left corner sensitivity mismatch: dNodesBdx(1,2:4,:) vs dNodesLdx(1,2:4,:).'); end
  if max(abs(sBR(:))) > tol, error('Bottom-right corner sensitivity mismatch: dNodesBdx(end,2:4,:) vs dNodesRdx(1,2:4,:).'); end
  if max(abs(sTL(:))) > tol, error('Top-left corner sensitivity mismatch: dNodesTdx(1,2:4,:) vs dNodesLdx(end,2:4,:).'); end
  if max(abs(sTR(:))) > tol, error('Top-right corner sensitivity mismatch: dNodesTdx(end,2:4,:) vs dNodesRdx(end,2:4,:).'); end
end

% PARAMETERS / GRIDS
u = linspace(0,1,Nu).';  v = linspace(0,1,Nv).';
[U,V] = ndgrid(u,v);                               % Nu x Nv

% COONS SURFACE (Nu * Nb * 3; 1 layer per component x,y,z)
t = reshape(t, Nu,1,3);
b = reshape(b, Nu,1,3);
l = reshape(l, 1,Nv,3);
r = reshape(r, 1,Nv,3);

b0 = b(1,:,:);  b1 = b(end,:,:);
t0 = t(1,:,:);  t1 = t(end,:,:);

R = (1-U).*l + U.*r + (1-V).*b + V.*t ...
  - (1-U).*(1-V).*b0 - (1-U).*V.*t0 - U.*(1-V).*b1 - U.*V.*t1;
  
% PACK NODES (u-first numbering)
ID = reshape(1:Nu*Nv,Nu,Nv);                       % IDs: u-first, then v
Nodes = [ID(:), reshape(R(:,:,1),[],1), reshape(R(:,:,2),[],1), reshape(R(:,:,3),[],1)];

% BOUNDARY NODES WITH NEW IDS
NodesB(:,1) = ID(:,1);                             % bottom: left -> right
NodesT(:,1) = ID(:,Nv);                            % top: left -> right
NodesL(:,1) = ID(1,:).';                           % left: bottom -> top
NodesR(:,1) = ID(Nu,:).';                          % right: bottom -> top

% SENSITIVITIES
dNodesdx = [];
if haveSens
  dNodesdx = zeros(Nu*Nv,4,nVar);

  dtdx = reshape(dNodesTdx(:,2:4,:), Nu,1,3,nVar);
  dbdx = reshape(dNodesBdx(:,2:4,:), Nu,1,3,nVar);
  dldx = reshape(dNodesLdx(:,2:4,:), 1,Nv,3,nVar);
  drdx = reshape(dNodesRdx(:,2:4,:), 1,Nv,3,nVar);

  db0dx = dbdx(1,:,:,:);
  db1dx = dbdx(end,:,:,:);
  dt0dx = dtdx(1,:,:,:);
  dt1dx = dtdx(end,:,:,:);

  dRdx = (1-U).*dldx + U.*drdx + (1-V).*dbdx + V.*dtdx ...
       - (1-U).*(1-V).*db0dx - (1-U).*V.*dt0dx - U.*(1-V).*db1dx - U.*V.*dt1dx;

  dNodesdx(:,2:4,:) = reshape(dRdx,Nu*Nv,3,nVar);
end
