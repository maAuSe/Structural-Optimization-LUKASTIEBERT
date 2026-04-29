function s = eltdef(Type)

switch lower(Type)
    case 'shell4'
        s.stress = stress_shell(4);
        
    case 'shell6'
        s.stress = stress_shell(3);
                
    case 'shell8'
        s.stress = stress_shell(4);
        
    case 'solid8'
        
        s.stress = stress_solid(8);
        
    case 'solid20'
        
        s.stress = stress_solid(8);
        
    case 'beam'
        
        s.stress = stress_bar;
    
    case 'truss'
        
        s.stress = stress_bar;
end


end

function stress = stress_shell(nnode)
stress.nstressc = 72;
stress.nnode = nnode; 
stress.sx.top = 6*(0:(nnode-1))+1;
stress.sx.mid = 6*(0:(nnode-1))+1+24;
stress.sx.bot = 6*(0:(nnode-1))+1+48;
stress.s3 = stress.sx;
stress.sy.top = 6*(0:(nnode-1))+1+1;
stress.sy.mid = 6*(0:(nnode-1))+1+1+24;
stress.sy.bot = 6*(0:(nnode-1))+1+1+48;
stress.s2 = stress.sy;
stress.sz.top = 6*(0:(nnode-1))+1+2;
stress.sz.mid = 6*(0:(nnode-1))+1+2+24;
stress.sz.bot = 6*(0:(nnode-1))+1+2+48;
stress.s1 = stress.sz;
stress.sxy.top = 6*(0:(nnode-1))+1+3;
stress.sxy.mid = 6*(0:(nnode-1))+1+3+24;
stress.sxy.bot = 6*(0:(nnode-1))+1+3+48;
stress.syz.top = 6*(0:(nnode-1))+1+4;
stress.syz.mid = 6*(0:(nnode-1))+1+4+24;
stress.syz.bot = 6*(0:(nnode-1))+1+4+48;
stress.sxz.top = 6*(0:(nnode-1))+1+5;
stress.sxz.mid = 6*(0:(nnode-1))+1+5+24;
stress.sxz.bot = 6*(0:(nnode-1))+1+5+48;
stress.header = [];
end

function stress = stress_solid(nnode)
stress.nstressc = 48;
stress.nnode = nnode;
stress.sx = (1:6:43);
stress.s3 = stress.sx;
stress.sy = (1:6:43)+1;
stress.s2 = stress.sy;
stress.sz = (1:6:43)+2;
stress.s1 = stress.sz;
stress.sxy = (1:6:43)+3;
stress.syz = (1:6:43)+4;
stress.sxz = (1:6:43)+5;
stress.header = [];
end

function stress = stress_bar
stress.nstressc = 1;
stress.sx = 1;
stress.s3 = stress.sx;
stress.sy = 1;
stress.s2 = stress.sy;
stress.sz = 1;
stress.s1 = stress.sz;
stress.syz = 1;
stress.sxy = 1;
stress.sxz = 1;
stress.header = [];
end