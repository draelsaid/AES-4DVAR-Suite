
function [x,xmat] = Mnl(x,fm,dt)  
     
    xmat(:,1) = x;
    
    xn = fm(xmat(:,1),dt);     
    xmat(:,2)  = x + xn;
     
    xn = fm(xmat(:,2),dt);
    xmat(:,3)  = x + xn/2;
     
    xn = fm(xmat(:,3),dt);
    xmat(:,4)  = x + xn/2;
     
    xn = fm(xmat(:,4),dt);
    xmat(:,5)  = x + xn;
     
    x  = (1/6)*( xmat(:,2) + 2*xmat(:,3) + 2*xmat(:,4) + xmat(:,5));

end