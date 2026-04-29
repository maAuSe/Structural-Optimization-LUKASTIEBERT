function F = tloads_beam(TLoad,Node,Section,Material)

%TLOADS_BEAM   Equivalent nodal forces for a beam element in the GCS.
%
%   F = tloads_beam(TLoad,Node,Section,Material)
%   computes the equivalent nodal forces of a temperature load 
%   (in the global coordinate system).
%
%   TLoad      Temperature loads        [dTm; dTy; dTz] (3 * 1)
%   Node       Node definitions         [x y z] (2 * 3)
%   Section    Section definitions      [A ky kz Ixx Iyy Izz yt yb zt zb]
%   Material   Material definitions     [E nu rho alpha]
%   F          Load vector              (12 * 1)
%
%   See also ELEMTLOADS, TLOADS_TRUSS.

% Pieter Reumers, Jef Wambacq
% 2018

A=Section(1);
Iyy=Section(5);
Izz=Section(6);
% Section heights in both y and z direction
if any(TLoad(2,:))
    try
        hy=Section(7)+Section(8);
    catch
        error('Properly define section dimensions yt and yb when applying a temperature gradient.');
    end
end
if any(TLoad(3,:))
    try
        hz=Section(9)+Section(10);
    catch
        error('Properly define section dimensions zt and zb when applying a temperature gradient.');
    end
end
E=Material(1);
if length(Material)>3
    alpha=Material(4);  % Material=[E nu rho alpha]
else 
    alpha=Material(3);  % Material=[E nu alpha]
end

dTm=TLoad(1,:);          % Average element temperature
dTy=TLoad(2,:);          % Temperature gradient in local y direction
dTz=TLoad(3,:);          % Temperature gradient in local z direction

Flcs=zeros(12,size(dTm,2));
% Axial forces
Flcs(1,:)=-E*A*alpha*dTm;
Flcs(7,:)= E*A*alpha*dTm;
% Bending moments Mz
if any(TLoad(2,:))
    Flcs(6,:)=  E*Izz*alpha*dTy/hy;
    Flcs(12,:)=-E*Izz*alpha*dTy/hy;
end
% Bending moments My
if any(TLoad(3,:))
    Flcs(5,:)= -E*Iyy*alpha*dTz/hz;
    Flcs(11,:)= E*Iyy*alpha*dTz/hz;
end

t=trans_beam(Node);
T=blkdiag(t,t,t,t);

F=T.'*Flcs; 
end
