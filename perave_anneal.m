% f=@(y) perave_opti(y(1)*1e-6,y(2),y(3),.66,2);
% simulannealbnd(f,[16,1.8,-20],[5,0,-80],[20,2*pi,-10]);

f=@(y) perave_opti(y(1)*1e-6,1.8,-20,y(2),2);
simulannealbnd(f,[16,.66],[5,0,.1],[20,2*pi,.99]);