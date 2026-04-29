function FLCS = loadslcs_shell2(DLoadLCS,L)

%LOADSLCS_SHELL2   Equivalent nodal forces for a SHELL2 element in the LCS.
%
%   FLCS = loadslcs_shell2(DLoadLCS,L)
%   computes the equivalent nodal forces of a distributed load
%   (in the local coordinate system).
%
%   DLoadLCS   Distributed loads        [n1localX; n1localY; n1localZ; ...]
%                                                                   (6 * 1)
%   L          Element length
%   FLCS       Load vector  (12 * 1)

% Mattias Schevenels
% April 2020

if size(DLoadLCS,1)==6 % only distributed forces
  px1 = DLoadLCS(1,:);
  py1 = DLoadLCS(2,:);
  mz1 = zeros(size(px1));
  px2 = DLoadLCS(4,:);
  py2 = DLoadLCS(5,:);
  mz2 = zeros(size(px1));
elseif size(DLoadLCS,1)==12 % distributed forces and moments
  px1 = DLoadLCS(1,:);
  py1 = DLoadLCS(2,:);
  mz1 = DLoadLCS(6,:);
  px2 = DLoadLCS(7,:);
  py2 = DLoadLCS(8,:);
  mz2 = DLoadLCS(12,:);
else
  error('SIZE(DloadLCS,1) must be 6 or 12.)');
end

if any(mz1~=mz2)
  error('SHELL2 element only supports constant distributed moments');
end
mz = mz1;

FLCS=zeros(6,size(DLoadLCS,2));

% loads in the local x-direction
FLCS(1,:)=L/6*(2*px1+px2);
FLCS(4,:)=L/6*(px1+2*px2);

% loads in the local y-direction
FLCS(2,:)=L/20*(7*py1+3*py2)-mz;
FLCS(5,:)=L/20*(3*py1+7*py2)+mz;
FLCS(3,:)=L^2/60*(3*py1+2*py2);
FLCS(6,:)=-L^2/60*(2*py1+3*py2);

