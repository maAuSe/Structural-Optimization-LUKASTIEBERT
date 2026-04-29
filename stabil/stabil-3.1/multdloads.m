function DLoads = multdloads(varargin)

%MULTDLOADS   Combine distributed loads.
%
%   DLoads = MULTDLOADS(DLoads_1,DLoads_2,...,DLoads_k)
%   combines the distributed loads of multiple load cases into one 3D array.
%   Each 3D-plane corresponds to a single load case. In the presence of 
%   partial distributed loads, every row will have identical starting 
%   and ending points in the 3rd dimension to allow for accurate combination 
%   of load cases. A distributed load with a starting and/or ending point 
%   value of 'NaN' will be considered as a load on the entire element
%   Calculations will be more efficient compared to a partial distributed 
%   load with starting point equal to zero and ending point equal to the 
%   length of the element.
%
%   DLoads_k   Distributed loads       [EltID n1globalX n1globalY n1globalZ ...]
%                                         without partial DLoads:  (nElem_k * 7)
%                                         with partial DLoads:     (nElem_k * 9)
%   DLoads     Distributed loads       [EltID n1globalX n1globalY n1globalZ ...]
%                                                             (maxnElem * 7 * k)
%                                                             (maxnElem * 9 * k)
%
%
%   See also ELEMLOADS.

% David Dooms, Wouter Dillen
% March 2009, December 2016

elemnumbers=[];

% get element indices with a distributed load
includePartialDLoads = 0;
for n=1:nargin
    DLoadn=varargin{n};
    elemnumbers=[elemnumbers; DLoadn(:,1)];
    if size(DLoadn,2)==9
        includePartialDLoads = 1;
    elseif size(DLoadn,2)~=7
        error('Input argument %i must have either 7 or 9 columns.',n);
    end
end
elemnumbers=unique(elemnumbers);


if includePartialDLoads == 0
    
    % Assemble DLoads
    DLoads=zeros(length(elemnumbers),7,nargin);
    DLoads(:,1,:)=repmat(elemnumbers,[1,1,nargin]);
    for n=1:nargin
        DLoadn=varargin{n};
        for i=1:size(DLoadn,1)
            irow = find(elemnumbers==DLoadn(i,1));
            DLoads(irow,2:7,n) = DLoads(irow,2:7,n) + DLoadn(i,2:7);
        end
    end
    
elseif includePartialDLoads == 1
    
    % Preallocate DLoads
    DLoadsonelement = zeros(length(elemnumbers),2,nargin);
    boolComplete    = cell(length(elemnumbers),nargin);
    complete        = cell(length(elemnumbers),nargin);
    Loc             = cell(length(elemnumbers),nargin);
    for n=1:nargin
        if size(varargin{n},2)==7
            varargin{n} = [varargin{n}, NaN(size(varargin{n},1),2)];
        end
        DLoadn = varargin{n};
        for i=1:length(elemnumbers)
            Loc{i,n} = find(DLoadn(:,1)==elemnumbers(i));
            boolComplete{i,n} = isnan(sum(DLoadn(Loc{i,n},8:9),2));
            complete{i,n} = (Loc{i,n}(boolComplete{i,n}));
            varargin{n}(complete{i,n},8:9) = NaN;
            DLoadsonelement(i,1,n) = length(complete{i,n});  % number of DLoads on the entire element
            DLoadsonelement(i,2,n) = length(Loc{i,n})-length(complete{i,n});  % number of partial DLoads on the element
        end
    end
    EltID = [];
    for i=1:length(elemnumbers)
        nrows = sum(DLoadsonelement(i,2,:),3);
        if sum(DLoadsonelement(i,1,:),3)>0
            nrows = nrows + 1;
        end
        EltID = [EltID; repmat(elemnumbers(i),[nrows,1])];
    end
    
    % Assemble DLoads:
    % * all DLoads on the entire element have a common global row entry between load cases
    % * all partial DLoads have unique entries to allow for accurate combination of load cases
    DLoads = zeros(length(EltID),9,nargin);
    DLoads(:,1,:) = repmat(EltID,[1,1,nargin]);
    for i=1:length(elemnumbers)  
        correction = 0;
        for n=1:nargin
            DLoadn = varargin{n};
            % collect DLoad information
            if isempty(complete{i,n})
                completeDLoad = [];
            else
                completeDLoad = sum(DLoadn(complete{i,n},2:end),1);
            end
            partial = Loc{i,n}(not(boolComplete{i,n}));
            partialDLoads = DLoadn(partial,2:end);
            % insert in DLoads
            irow = find(elemnumbers(i)==DLoads(:,1),1,'first');
            if sum(DLoadsonelement(i,1,:),3)>0
                if ~isempty(completeDLoad)
                    DLoads(irow,2:end,n) = completeDLoad;
                    DLoads(irow,8:9,:) = NaN;
                end
                irow = irow + 1;
            end
            if ~isempty(partialDLoads)
                irow = irow + correction;
                DLoads(irow:irow+size(partialDLoads,1)-1,2:end,n) = partialDLoads;
                correction = correction + size(partialDLoads,1);
            end
        end
    end
end

