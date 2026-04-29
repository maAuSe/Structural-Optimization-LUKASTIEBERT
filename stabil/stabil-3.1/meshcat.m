function [Nodes,Elements] = meshcat(varargin)

if nargin>4
  [Nodes1,Elements1] = meshcat(varargin{1:4});
  [Nodes,Elements] = meshcat(Nodes1,Elements1,varargin{5:end});
elseif nargin==4
  Nodes1 = varargin{1};
  Elements1 = varargin{2};
  Nodes2 = varargin{3};
  Elements2 = varargin{4};

  N = max(Nodes1(:,1));
  Nodes2(:,1) = Nodes2(:,1)+N;
  Elements2(:,5:end) = Elements2(:,5:end)+N;
  Nodes = [Nodes1; Nodes2];
  Elements = [Elements1; Elements2];

  for k = numel(Elements(:,1:4))+1:numel(Elements)
    Elements(k) = find(Nodes(:,1)==Elements(k));
  end
  Nodes(:,1) = [1:size(Nodes,1)]';

  [~,ia,ic] = uniquetol(Nodes(:,2:4)+2*max(max(abs(Nodes(:,2:4)))),1e-8,'ByRows',true);

  Nodes = Nodes(ia,:);
  Elements(:,5:end) = ic(Elements(:,5:end));

  Nodes(:,1) = [1:size(Nodes,1)]';

  Elements(:,1) = [1:size(Elements,1)]';

else
  error('Incorrect number of input arguments.');
end
