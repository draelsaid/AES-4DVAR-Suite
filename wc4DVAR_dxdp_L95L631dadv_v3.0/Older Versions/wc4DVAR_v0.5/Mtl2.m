
function [x,dx] = Mtl2(dx,x,fm,dt)    
    dxn = fm(x,dx);     
     x  = x + dt*xn;
     dx = dx + dt*dxn;
end