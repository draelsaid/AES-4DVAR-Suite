function dx = Madjrk2(dx,x,fadj,fnl,dt)  
  [~,xmat] = Mnl(x,fnl,dt);   % Call to non-linear model, to import linearisation states (k1,k2,k3,k4)
  
 zz = fadj(xmat(:,1),dx,dt);
 dx = dx + zz;
 
end