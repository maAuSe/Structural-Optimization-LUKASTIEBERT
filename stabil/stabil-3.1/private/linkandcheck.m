function p3=linkandcheck(p1,p2,name)

%LINKANDCHECK   Link two matrices and check for presence.
%
%   p3=linkandcheck(p1,p2,name)
%   Link two matrices and check for presence
%
%   p1     Vector with elements to be link and checked
%   p2     Vector with elements to be link and checked
%   name   Name of elements in p2
%   p3     Lowest absolute index in p2 for each element in p1

[~,p3] = ismember(round(p1),round(p2));
loc = find(p3==0);
if ~isempty(loc)
    error('%s %i is not defined.',name,p1(loc(1)))
end
end
