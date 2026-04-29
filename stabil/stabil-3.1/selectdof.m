function [L,I]=selectdof(DOF,seldof,varargin)

%SELECTDOF   Select degrees of freedom.
%
%       L=selectdof(DOF,seldof)
%   [L,I]=selectdof(DOF,seldof)
%       L=selectdof(DOF,seldof,'Ordering',ordering)
%   creates the matrix to extract degrees of freedom from the global degrees of 
%   freedom by matrix multiplication.
%
%   DOF        Degrees of freedom  (nDOF * 1)
%   seldof     Selected DOF labels (kDOF * 1)
%   L          Selection matrix (kDOF * nDOF)
%   I          Index vector (kDOF * 1)
%   ordering   'seldof','DOF' or 'sorted'                      Default: 'seldof'
%              Ordering of L and I similar as seldof, DOF or sorted  

%   See also UNSELECTDOF.

% David Dooms (March 2008), Miche Jansen (January 2012)


% PREPROCESSING
DOF=round(DOF(:)*100);
seldof=round(seldof(:)*100);
if nargin<3
    varargin={};
end
paramlist=varargin;
[ordering,paramlist]=cutparam('Ordering','seldof',paramlist,{'seldof','DOF','sorted'});
if ~isempty(paramlist)
    error(['Undefined parameter ''' paramlist{1} '''.']); end
if ~ isempty(find(seldof==0.00,1))
    error('The wild card 0.00 is not allowed')
end

nDOF=length(DOF);
kDOF=length(seldof);

% DOF should be a unique vector
checkunique(DOF,'Degree of freedom');
% [DOF,dum,jsel] = unique(DOF);

% expand seldof
% check wild cards
wildcards = zeros(kDOF,1);
wildcards(rem(seldof,100)==0) = 1;
wildcards(floor(seldof/100)==0) = 2;
if any(wildcards>0)    
[useldof,dum,isel] = unique(seldof);
[al,a]=ismember(DOF,useldof);
[bl,b]=ismember(floor(DOF/100)*100,useldof);
[cl,c]=ismember(rem(DOF,100),useldof);
seldofx=[DOF(al);DOF(bl);DOF(cl)]; %expanded seldof
[dum,ind]=sort([a(al);b(bl);c(cl)]);
ind = mat2cell(ind,accumarray(dum,1,[length(useldof),1]));
seldofx2 = cell(length(ind),1);
for it = 1:length(ind)
   if isempty(ind{it})
     seldofx2{it} = useldof(it); 
   else
     seldofx2{it} = seldofx(ind{it});  
   end
end
seldofx = cell2mat(seldofx2(isel));
else
seldofx = seldof;    
end


[indl,indj] = ismember(seldofx,DOF);
switch lower(ordering)
    case 'seldof'
% I = indj(indl);
    I = indj;
    indi=(1:length(indj)).';
    kDOF = length(indj);
    indi = indi(indl);
    indj = indj(indl);
    case 'dof'
    indj = unique(indj(indl));
    I = indj;
    kDOF = length(indj);
    indi=(1:kDOF).';    
    case 'sorted'
     indj = indj(indl);
    [sortedDOF,indd] = unique(DOF(indj),'first');
    indi=(1:length(indj)).';
    indj = indj(indd);
    indi = indi(indd);
    I = indj;
    kDOF = length(indj);
%     jsel = sort(jsel);
end

L = sparse(indi,indj,1,kDOF,nDOF);
% L = L(:,jsel);
end

%-------------------------------------------------------------------------------
% CUT PARAMETER FROM LIST
function [value,paramlist]=cutparam(name,default,paramlist,allowed)
value=default;
for iarg=length(paramlist)-1:-1:1
  if strcmpi(name,paramlist{iarg})
    value=paramlist{iarg+1};
    paramlist=paramlist([1:iarg-1 iarg+2:end]);
    break
  end
end
if nargin==4
  if isempty(find(strcmpi(value,allowed)))
    error(['Parameter ' name ' should be one of ' sprintf('''%s'' ',allowed{:}) '.']);
  end
end
end