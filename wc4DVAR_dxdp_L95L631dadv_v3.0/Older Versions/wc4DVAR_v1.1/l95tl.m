function [dy] = l95tl(x,dx,dt)

 N = length(x);     % Number of variables automatically set by user
 
 for i=1:N
 % Periodic (cyclic) Boundary Conditions    
    ip1 = i+1; if ip1>N, ip1 = ip1-N; end 
    im1 = i-1; if im1<1, im1 = im1+N; end
    im2 = i-2; if im2<1, im2 = im2+N; end
 % RHS of L95 TL
   dy(i,1) = dt*(- dx(im2)*x(im1) + dx(im1)*(x(ip1) - x(im2)) - dx(i) + dx(ip1)*x(im1)); 
 end
 
end