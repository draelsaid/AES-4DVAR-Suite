
function [x] = Mnl2(x,fm,dt)  
     
   xn = fm(x);     
   x  = x + dt*xn;

end