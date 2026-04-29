function [value,paramlist]=cutparam(name,default,paramlist);
% CUT PARAMETER FROM LIST
value=default;
for iarg=length(paramlist)-1:-1:1
    if strcmpi(name,paramlist{iarg})
        value=paramlist{iarg+1};
        paramlist=paramlist([1:iarg-1 iarg+2:end]);
        break
    end
end