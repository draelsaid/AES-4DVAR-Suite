function [dy,M] = l95tlmtx(x,dx,dt)

N = length(x);  % Number of variables automatically set by user

% L95 Matrix (for cond number calc)
M = zeros(N);
M(1,1) = x(1); M(1,2) = x(N);
M(2,1) = x(3)-x(N); M(2,2) = x(2); M(2,3) = x(1); 
for j =3:N-1
 M(j,j-2) = -x(j-1);
 M(j,j-1) = x(j+1) - x(j-2);
 M(j,j) = x(j);
 M(j,j+1) = x(j-1);
end
M(N,N-2) = x(N-1); M(N,N-1) = x(1) - x(N-2); M(N,N) = x(1);
 
 for i=1:N
 % Periodic (cyclic) Boundary Conditions    
    ip1 = i+1; if ip1>N, ip1 = ip1-N; end 
    im1 = i-1; if im1<1, im1 = im1+N; end
    im2 = i-2; if im2<1, im2 = im2+N; end
 % RHS of L95 TL
   dy(i,1) = dt*(- dx(im2)*x(im1) + dx(im1)*(x(ip1) - x(im2)) - dx(i) + dx(ip1)*x(im1)); 
 end
 
end