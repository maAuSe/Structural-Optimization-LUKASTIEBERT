function checkunique(p,name)

%CHECKUNIQUE   Check if vector contains unique elements.
%
%   checkunique(p,name)
%   find non-unique elements in a vector and print error
%
%   p      Vector with elements to be checked
%   name   Name of elements in p
 
if length(p)~=length(unique(p))
   for ip=1:length(p)
      loc=find(p==p(ip));
      if length(loc)>1
         error('%s %g is multiply defined.',name,p(ip))
      end
   end
end
end
