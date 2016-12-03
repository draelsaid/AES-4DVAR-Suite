
function [x,xmat] = Mnlrk2(x,fm,dt)  
     
    xmat(:,1) = x;
    
    xn = fm(x,dt);     
    x  = x + xn;

end