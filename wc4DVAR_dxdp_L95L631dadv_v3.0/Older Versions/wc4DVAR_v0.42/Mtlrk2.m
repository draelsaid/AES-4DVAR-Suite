function [dx] = Mtlrk2(dx,x,ftl,fnl,dt)  
   
  [~,xmat] = Mnl(x,fnl,dt);  
  
   [dxn] = ftl(xmat(:,1),dx,dt);      
     dx = dx + dxn;
    
end
