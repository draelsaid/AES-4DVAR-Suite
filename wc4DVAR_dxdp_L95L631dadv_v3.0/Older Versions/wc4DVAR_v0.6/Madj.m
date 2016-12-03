function dx = Madj(dx,x,fadj,fnl,dt)  
  [~,xmat] = Mnl(x,fnl,dt);   % Call to non-linear model, to import linearisation states (k1,k2,k3,k4)
  
 xx = dx/6;
 qq = xx;

 zz = fadj(xmat(:,4),qq,dt);
 dx = dx + zz;
 qq = zz;
  
 qq = qq + 2*xx;
 zz = fadj(xmat(:,3),qq,dt);
 dx = dx + zz;
 qq = 0.5*zz;
     
 qq = qq + 2*xx;
 zz = fadj(xmat(:,2),qq,dt);
 dx = dx + zz;
 qq = 0.5*zz;
     
 qq = qq + xx;
 zz = fadj(xmat(:,1),qq,dt);
 dx = dx + zz;
 
end