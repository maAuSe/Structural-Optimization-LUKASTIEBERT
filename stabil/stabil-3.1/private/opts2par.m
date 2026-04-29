function par=opts2par(Options,field,default)
%% Process options
if isfield(Options,field)
  par=getfield(Options,field);
else
  par=default;
end
end