function y = l95(x,dt)

N = length(x); % Number of variables automatically set by user
F=8;           % Forcing term = 8 for chaos

 for i=1:N 
% Periodic (cyclic) Boundary Conditions    
    ip1 = i+1; if ip1>N, ip1 = ip1-N; end 
    im1 = i-1; if im1<1, im1 = im1+N; end
    im2 = i-2; if im2<1, im2 = im2+N; end
% RHS of L95
 y(i,1) = dt*(-x(im2)*x(im1) + x(im1)*x(ip1) - x(i) + F); 
 end

end