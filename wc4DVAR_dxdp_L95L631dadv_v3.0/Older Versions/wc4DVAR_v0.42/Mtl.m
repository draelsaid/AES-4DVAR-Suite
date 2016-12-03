function [dx] = Mtl(dx,x,ftl,fnl,dt)  
   
  [~,xmat] = Mnl(x,fnl,dt);  
  
   [dxn] = ftl(xmat(:,1),dx,dt);      
     dk1 = dx + dxn;
     
   [dxn] = ftl(xmat(:,2),dk1,dt);
     dk2 = dx + dxn/2;
     
   [dxn] = ftl(xmat(:,3),dk2,dt);
     dk3 = dx + dxn/2;
     
   [dxn] = ftl(xmat(:,4),dk3,dt);
     dk4 = dx + dxn;
     
 dx = (dk1+2*dk2+2*dk3+dk4)/6;

end
