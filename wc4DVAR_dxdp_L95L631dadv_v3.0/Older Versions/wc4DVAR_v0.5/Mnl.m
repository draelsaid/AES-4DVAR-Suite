
function [x] = Mnl(x,fm,dt)  
     
    [xn] = fm(x);     
     k1  = x + dt*xn;
     
    [xn] = fm(k1);
     k2  = x + (dt/2)*xn;
     
    [xn] = fm(k2);
     k3  = x + (dt/2)*xn;
     
    [xn] = fm(k3);
     k4  = x + dt*xn;
     
    x  = (1/6)*(k1+2*k2+2*k3+k4);

end