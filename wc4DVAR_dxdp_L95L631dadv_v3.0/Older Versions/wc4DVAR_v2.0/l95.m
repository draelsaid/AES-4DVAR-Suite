function y = l95(x,dt)

N = length(x); % Number of variables automatically set by user
F=8;           % Forcing term = 8 for chaos

% L95 Matrix (for cond number calc)
M(1,1) = x(1); M(1,2) = x(N);
M(2,1) = x(3)-x(N); M(2,2) = x(2); M(2,3) = x(1); 
for j =3:N-1
 M(j-2,j) = -x(j-1);
 M(j-1,j) = x(j+1) - x(j-2);
 M(j,j) = x(j);
 M(j+1,j) = x(j-1);
end
M(N,N-2) = x(N-1); M(N,N-1) = x(1) - x(N-2); M(N,N) = x(1);

 for i=1:N 
% Periodic (cyclic) Boundary Conditions    
    ip1 = i+1; if ip1>N, ip1 = ip1-N; end 
    im1 = i-1; if im1<1, im1 = im1+N; end
    im2 = i-2; if im2<1, im2 = im2+N; end
% RHS of L95
 y(i,1) = dt*(-x(im2)*x(im1) + x(im1)*x(ip1) - x(i) + F); 
 end
