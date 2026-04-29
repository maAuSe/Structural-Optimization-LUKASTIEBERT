function [DLoadLCS,dDLoadLCSdx] = dloadgcs2lcs(T,DLoad,dTdx,dDLoaddx)

%DLOADGCS2LCS   Distributed load transformation from GCS to LCS.
%
%   DLoadLCS = DLOADGCS2LCS(T,DLoad)
%       transforms the distributed load definitions in the global
%       coordinate system (algebraic convention) to the local coordinate
%       system (beam convention).
%
%   [DLoadLCS,dDLoadLCSdx] = DLOADGCS2LCS(T,DLoad,dTdx,dDLoaddx)
%       transforms the distributed load definitions in the global
%       coordinate system (algebraic convention) to the local coordinate
%       system (beam convention), and additionally computes the derivatives
%       of the distributed load information with respect to the design
%       variables x.
%
%   T           Element transformation matrix        (6 * 6)
%   DLoad       Distributed loads in GCS             [n1globalX; n1globalY; n1globalZ; ...]  
%               	(6/8 * nLC * nDLoads)
%   dTdx        Transformation matrix derivatives    (6 * 6 * nVar)
%   dDLoaddx    Distributed loads derivatives (GCS)  (SIZE(DLoad) * nVar)
%   DLoadLCS    Distributed loads in LCS             [n1localX; n1localY; n1localZ; ...]  
%               	(6/8 * nLC * nDLoads)
%   dDLoadLCSdx Distributed loads derivatives (LCS)  (SIZE(DLoadLCS) * nVar)    

% Wouter Dillen
% February 2017

if nargin<3, dTdx = []; end
if nargin<4, dDLoaddx = []; end

nVar = 0;
if nargout>1 && (~isempty(dTdx) || ~isempty(dDLoaddx))
    nVar = max(size(dTdx,3),size(dDLoaddx,3));
end

if nVar==0 || isempty(dTdx), dTdx = zeros([size(T),nVar]); end
if nVar==0 || isempty(dDLoaddx), dDLoaddx = zeros(size(DLoad,1),size(DLoad,2),nVar,size(DLoad,3)); end


nDLoadsOnElement=size(DLoad,3);

% DLoadLCS
DLoadLCS = DLoad;
for i=1:nDLoadsOnElement
    DLoadLCS(1:6,:,i) = T*DLoad(1:6,:,i);
end
   
% dDLoadLCSdx
dDLoadLCSdx = zeros(size(DLoad,1),size(DLoad,2),nVar,size(DLoad,3));
for n=1:nVar
    
    %any(reshape(dTdx(:,:,n),[],1)~=0)
    %any(reshape(dDLoaddx(1:6,:,n,:),[],1)~=0)
    
    % dTdx(:,:,n) typisch maar bij nVar niet volledig 0
    % dDLoaddx(1:6,:,:,:)   bijna altijd gelijk aan 0
    % dDLoaddx(7:8,:,:,:)   ook dikwijls 0 of NaN, tenzij
    %                       bij 2 ontwerpvariabelen (a of b).
    %       maakt dat eigenlijk nog wel uit als dDLoaddx al gelijk is aan 0?
    
    
    % be aware of numerical inaccuracies in dDLoadLCSdx when a or b
    % approaches (but is not equal to) zero

    if any(reshape(dTdx(:,:,n),[],1)~=0)
        
        for i=1:nDLoadsOnElement         
            dDLoadLCSdx(1:6,:,n,i)=dTdx(:,:,n)*DLoad(1:6,:,i)+T*dDLoaddx(1:6,:,n,i);
        end
        if size(DLoad,1)==8  
            % unit vector along the local x-axis (Nx)
            t=T(1,1:3); dtdx=dTdx(1,1:3,n);
            for i=1:nDLoadsOnElement         
                % a
                a = DLoad(7,:,i);
                if all(isnan(a))
                    dDLoadLCSdx(7,:,n,i) = NaN;
                else
                    dadx = dDLoaddx(7,:,n,i);
                    loc = find(dadx~=0); 
                    if ~isempty(loc)
                        dDLoadLCSdx(7,loc,n,i) = diag((dadx(loc).'*t+a(loc).'*dtdx)*(t.'*a(loc))).'./a(loc);
                    end
                end 
                % b
                b = DLoad(8,:,i);
                if all(isnan(b))
                    dDLoadLCSdx(8,:,n,i) = NaN;
                else
                    dbdx = dDLoaddx(8,:,n,i);
                    loc = find(dbdx~=0); 
                    if ~isempty(loc)
                        dDLoadLCSdx(8,loc,n,i) = diag((dbdx(loc).'*t+b(loc).'*dtdx)*(t.'*b(loc))).'./b(loc);
                    end
                end 
            end
        end
        
    else
        
        for i=1:nDLoadsOnElement         
            dDLoadLCSdx(1:6,:,n,i)=T*dDLoaddx(1:6,:,n,i);
        end
        if size(DLoad,1)==8  
            % unit vector along the local x-axis (Nx)
            t=T(1,1:3);
            for i=1:nDLoadsOnElement         
                % a
                a = DLoad(7,:,i);
                if all(isnan(a))
                    dDLoadLCSdx(7,:,n,i) = NaN;
                else
                    dadx = dDLoaddx(7,:,n,i);
                    loc = find(dadx~=0); 
                    if ~isempty(loc)
                        dDLoadLCSdx(7,loc,n,i) = diag((dadx(loc).'*t)*(t.'*a(loc))).'./a(loc);
                    end
                end 
                % b
                b = DLoad(8,:,i);
                if all(isnan(b))
                    dDLoadLCSdx(8,:,n,i) = NaN;
                else
                    dbdx = dDLoaddx(8,:,n,i);
                    loc = find(dbdx~=0); 
                    if ~isempty(loc)
                        dDLoadLCSdx(8,loc,n,i) = diag((dbdx(loc).'*t)*(t.'*b(loc))).'./b(loc);
                    end
                end 
            end
        end
        
    end
end


