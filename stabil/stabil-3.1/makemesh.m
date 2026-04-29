function [Nodes,Elements,Edge1,Edge2,Edge3,Edge4,Normals] = makemesh(Line1,Line2,Line3,Line4,m,n,Type,Section,Material,varargin)

%MAKEMESH   Creates a mesh of quadrilateral elements for a surface defined by 4 lines
%
%   [Nodes,Elements,Edge1,Edge2,Edge3,Edge4,Normals] =
%   makemesh(Line1,Line2,Line3,Line4,m,n,Type,Section,Material)
%
%   Line1,2,3,4 Lines that define the surface. (n1 * 3)
%               These lines:
%                   -should be defined in a clockwise or counter-clockwise order
%                    (this means: first point of Line2 is the last point of Line1
%                    in the case of a counter-clockwise direction)
%                   -can be defined by any number of points
%
%   m           Number of elements that divide Line2 and Line4
%   n           Number of elements that divide Line1 and Line3
%   Type        Cell containing Type number and name        {TypeID EltName}
%   Section     Section ID
%   Material    Material ID
%
%   Nodes       Nodes definitions       [NodID x y z]
%   Elements    Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Edge1,2,3,4 vector containing Nodes on Line1,2,..  ((m+1) * 1) or ((n+1) * 1)
%   Normals     Normals of shell surface in Nodes [NodID nx ny nz]  
%
%   makemesh(...,ParamName,ParamValue)sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'L1method'   Specify the interpolation scheme for Line1.
%                default: 'spline'. This method is used in the MATLAB
%                INTERP1 function. See DOC INTERP1 for available
%                interpolation schemes.
%   'L2method'
%   'L3method'
%   'L4method'
%   'UniformSpacing'    {'on',default: 'off'}
%               If UniformSpacing is set to 'off', the distance between
%               consecutive nodes along the edges of the mesh is based on the spacing
%               between the points in the corresponding Line.
%               If UniformSpacing is 'on', makemesh attempts to distribute
%               the nodes uniformly along the lines regardless of the
%               initial spacing between the points by approximating the
%               natural parametrizations of the lines.
%
%
%   See also GRID_SHELL8, GRID_SHELL4.

% Miche Jansen
% 2009

% PARAMETERS

if nargin > 9
    paramlist=varargin;
else paramlist = [];
end

if nargout > 6; lnormals = true; else lnormals = false; end

[LinterpMode{1},paramlist]=cutparam('L1method','spline',paramlist);
[LinterpMode{2},paramlist]=cutparam('L2method','spline',paramlist);
[LinterpMode{3},paramlist]=cutparam('L3method','spline',paramlist);
[LinterpMode{4},paramlist]=cutparam('L4method','spline',paramlist);
[UniformSpacing,paramlist]=cutparam('UniformSpacing','off',paramlist);


TypeName = Type{1,2};

[s,t,NodeNum,Elements] = eval(['grid_' TypeName '(m,n,Type,Section,Material)']);

Lines = {Line1;Line2;Line3;Line4};
% sn = {s;t;-s;-t};
sn = {s;t;s;t};
spdir = [1;1;-1;-1]; 
L = cell(4,1);

for iLine=1:4
    nsp = size(Lines{iLine},1);
    sp = spdir(iLine)*linspace(-1,1,nsp);
    
    if strcmpi(UniformSpacing,'on')
        sp2 = spdir(iLine)*linspace(-1,1,nsp*10);
        Lines{iLine} = interp1(sp,Lines{iLine},sp2,LinterpMode{iLine});
        dL2 =  diff(Lines{iLine}); dL2 = sqrt(sum(dL2.^2,2));
        sp = 2*[0;cumsum(dL2)]/sum(dL2)-1;
        sp = spdir(iLine)*sp;
    end
    
    L{iLine} = interp1(sp,Lines{iLine},sn{iLine},LinterpMode{iLine}).';
    
    if lnormals
       dL{iLine}  = bsxfun(@rdivide,diff(L{iLine},1,2),diff(sn{iLine}));
       dL{iLine}  = conv2(dL{iLine},[1,1]/2); % extrapolate to nodes
    end
end

N1_s = (1-s)/2;
N1_t = (1-t)/2;
N2_s = (1+s)/2;
N2_t = (1+t)/2;

P_tx = N1_t.'*L{1}(1,:)+N2_t.'*L{3}(1,:);
P_ty = N1_t.'*L{1}(2,:)+N2_t.'*L{3}(2,:);
P_tz = N1_t.'*L{1}(3,:)+N2_t.'*L{3}(3,:);

P_sx = L{4}(1,:).'*N1_s+L{2}(1,:).'*N2_s;
P_sy = L{4}(2,:).'*N1_s+L{2}(2,:).'*N2_s;
P_sz = L{4}(3,:).'*N1_s+L{2}(3,:).'*N2_s;

% de volgende formules gaan ervan uit dat de beginpunten van lijn 1 en 3
% gelijk zijn aan de eindpunten van lijn 2 en 4 en omgekeerd
P_stx = N1_t.'*N1_s*Line1(1,1)+N2_t.'*N1_s*Line3(end,1)+ ...
    N1_t.'*N2_s*Line1(end,1)+N2_t.'*N2_s*Line3(1,1);
P_sty = N1_t.'*N1_s*Line1(1,2)+N2_t.'*N1_s*Line3(end,2)+ ...
    N1_t.'*N2_s*Line1(end,2)+N2_t.'*N2_s*Line3(1,2);
P_stz = N1_t.'*N1_s*Line1(1,3)+N2_t.'*N1_s*Line3(end,3)+ ...
    N1_t.'*N2_s*Line1(end,3)+N2_t.'*N2_s*Line3(1,3);

phi_x = P_tx+P_sx-P_stx;
phi_y = P_ty+P_sy-P_sty;
phi_z = P_tz+P_sz-P_stz;

fnod = ~isnan(NodeNum);
Nodes = [NodeNum(fnod(:)),phi_x(fnod(:)),phi_y(fnod(:)),phi_z(fnod(:))];

if nargout > 2
    Edge4 = NodeNum(:,1);
    Edge1 = NodeNum(end,:).';
    Edge2 = NodeNum(:,end);
    Edge3 = NodeNum(1,:).';
end


% Compute normals on the surface in the nodes
if lnormals
    dN1_ds = -ones(size(s))/2;
    dN1_dt = -ones(size(t))/2;
    dN2_ds = ones(size(s))/2;
    dN2_dt = ones(size(t))/2;
    
    dP_tx_dt = dN1_dt.'*L{1}(1,:)+dN2_dt.'*L{3}(1,:);
    dP_ty_dt = dN1_dt.'*L{1}(2,:)+dN2_dt.'*L{3}(2,:);
    dP_tz_dt = dN1_dt.'*L{1}(3,:)+dN2_dt.'*L{3}(3,:);
    dP_tx_ds = N1_t.'*dL{1}(1,:)+N2_t.'*dL{3}(1,:);
    dP_ty_ds = N1_t.'*dL{1}(2,:)+N2_t.'*dL{3}(2,:);    
    dP_tz_ds = N1_t.'*dL{1}(3,:)+N2_t.'*dL{3}(3,:);

    dP_sx_ds = L{4}(1,:).'*dN1_ds+L{2}(1,:).'*dN2_ds;
    dP_sy_ds = L{4}(2,:).'*dN1_ds+L{2}(2,:).'*dN2_ds;    
    dP_sz_ds = L{4}(3,:).'*dN1_ds+L{2}(3,:).'*dN2_ds;
    dP_sx_dt = dL{4}(1,:).'*N1_s+dL{2}(1,:).'*N2_s;
    dP_sy_dt = dL{4}(2,:).'*N1_s+dL{2}(2,:).'*N2_s;
    dP_sz_dt = dL{4}(3,:).'*N1_s+dL{2}(3,:).'*N2_s;
    
    % de volgende formules gaan ervan uit dat de beginpunten van lijn 1 en 3
    % gelijk zijn aan de eindpunten van lijn 2 en 4 en omgekeerd
    dP_stx_ds = N1_t.'*dN1_ds*Line1(1,1)+N2_t.'*dN1_ds*Line3(end,1)+ ...                                     
        N1_t.'*dN2_ds*Line1(end,1)+N2_t.'*dN2_ds*Line3(1,1);                                             
    dP_sty_ds = N1_t.'*dN1_ds*Line1(1,2)+N2_t.'*dN1_ds*Line3(end,2)+ ...                                     
        N1_t.'*dN2_ds*Line1(end,2)+N2_t.'*dN2_ds*Line3(1,2);                                             
    dP_stz_ds = N1_t.'*dN1_ds*Line1(1,3)+N2_t.'*dN1_ds*Line3(end,3)+ ...                                     
        N1_t.'*dN2_ds*Line1(end,3)+N2_t.'*dN2_ds*Line3(1,3);           
    dP_stx_dt = dN1_dt.'*N1_s*Line1(1,1)+dN2_dt.'*N1_s*Line3(end,1)+ ...                                     
        dN1_dt.'*N2_s*Line1(end,1)+dN2_dt.'*N2_s*Line3(1,1);                                             
    dP_sty_dt = dN1_dt.'*N1_s*Line1(1,2)+dN2_dt.'*N1_s*Line3(end,2)+ ...                                     
        dN1_dt.'*N2_s*Line1(end,2)+dN2_dt.'*N2_s*Line3(1,2);                                             
    dP_stz_dt = dN1_dt.'*N1_s*Line1(1,3)+dN2_dt.'*N1_s*Line3(end,3)+ ...                                     
        dN1_dt.'*N2_s*Line1(end,3)+dN2_dt.'*N2_s*Line3(1,3);                                             
    
    dphi_x_ds = dP_tx_ds+dP_sx_ds-dP_stx_ds;
    dphi_y_ds = dP_ty_ds+dP_sy_ds-dP_sty_ds;
    dphi_z_ds = dP_tz_ds+dP_sz_ds-dP_stz_ds;
    dphi_x_dt = dP_tx_dt+dP_sx_dt-dP_stx_dt;
    dphi_y_dt = dP_ty_dt+dP_sy_dt-dP_sty_dt;
    dphi_z_dt = dP_tz_dt+dP_sz_dt-dP_stz_dt;
    % normal -> cross product 
    normal_x = dphi_y_ds.*dphi_z_dt-dphi_y_dt.*dphi_z_ds;
    normal_y = -(dphi_x_ds.*dphi_z_dt-dphi_z_ds.*dphi_x_dt);
    normal_z = dphi_x_ds.*dphi_y_dt-dphi_y_ds.*dphi_x_dt;
    nnormal = sqrt(normal_x.^2+normal_y.^2+normal_z.^2);
    normal_x = normal_x./nnormal;
    normal_y = normal_y./nnormal;
    normal_z = normal_z./nnormal;    
    Normals = [NodeNum(fnod(:)),normal_x(fnod(:)),normal_y(fnod(:)),normal_z(fnod(:))];   
end