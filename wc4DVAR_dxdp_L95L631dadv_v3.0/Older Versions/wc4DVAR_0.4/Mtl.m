
function [x,dx] = Mtl(dx,x,fm,dt)  
     
    [xn,dxn] = fm(x,dx);     
     k1  = x + dt*xn;
     dk1 = dx + dt*dxn;
     
    [xn,dxn] = fm(k1,dk1);
     k2  = x + (dt/2)*xn;
     dk2 = dx + (dt/2)*dxn;
     
    [xn,dxn] = fm(k2,dk2);
     k3  = x + (dt/2)*xn;
     dk3 = dx + (dt/2)*dxn;
     
    [xn,dxn] = fm(k3,dk3);
     k4  = x + dt*xn;
     dk4 = dx + dt*dxn;
     
    x  = (1/6)*(k1+2*k2+2*k3+k4);
     dx = (1/6)*(dk1+2*dk2+2*dk3+dk4);

end
