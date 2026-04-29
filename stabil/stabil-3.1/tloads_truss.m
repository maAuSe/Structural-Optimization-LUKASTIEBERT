function F = tloads_truss(TLoad,Node,Section,Material)

%TLOADS_TRUSS   Equivalent nodal forces for a truss element in the GCS.
%
%   F = tloads_truss(TLoad,Node,Section,Material)
%   computes the equivalent nodal forces of a temperature load 
%   (in the global coordinate system).
%
%   TLoad      Temperature loads        [n1globalX; n1globalY; n1globalZ; ...]
%                                                                   (6 * 1)
%   Node       Node definition          [x y z] (2 * 3)
%   Section    Section definition       [A ...]
%   Material   Material definition      [E nu rho alpha]
%   F          Load vector              (6 * 1)
%
%   See also ELEMTLOADS, TLOADS_BEAM.

% Jef Wambacq, Pieter Reumers
% 2018

if any(TLoad(2,:)) || any(TLoad(3,:))
    warning('Temperature gradient for truss elements is not possible. Using average temperature only.');
end

A=Section(1);
E=Material(1);
if length(Material)>3
    alpha=Material(4);  % Material=[E nu rho alpha]
else 
    alpha=Material(3);  % Material=[E nu alpha]
end

Tm=TLoad(1,:);           % Average element temperature

Flcs=zeros(6,size(Tm,2));
Flcs(1,:)=-E*A*alpha*Tm;
Flcs(4,:)= E*A*alpha*Tm;

t=trans_truss(Node);
T=blkdiag(t,t);

F=T.'*Flcs; 
end
