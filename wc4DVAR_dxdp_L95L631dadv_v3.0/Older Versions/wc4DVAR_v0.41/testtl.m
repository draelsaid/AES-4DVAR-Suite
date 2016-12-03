% Testing the tangent-linear model with the numerical integration scheme
clear all
close all

N=40; dt=0.025; n=200;

fnl = @l95;
schnl = @rk4nl;

ftl = @l95tl;
schtl = @rk4tl;

% Correctness Test
for i=1:16
alpha = 10^(1-i);
dx = alpha*randn(N,1);
x = randn(N,1);
xdx = x + dx;

Mxdx = Mtl(dx,x,ftl,fnl,schtl,schnl,dt,n); 
mxdx = Mnl(xdx,fnl,schnl,dt,n); mxdx = mxdx(:,n+1);
mx   = Mnl(x,fnl,schnl,dt,n); mx = mx(:,n+1);

alphadx(i)  = norm(mxdx - mx)/norm(Mxdx);          % ~1
alphadx2(i) = alphadx(i)-1;                        % ~10^(-8)
alphadx3(i) = norm(mxdx - mx - Mxdx)/norm(Mxdx);   % ~10^(-12)
end

% Validity Test
dx = alpha*randn(N,1);
x = randn(N,1);
xdx = x + dx;

Mxdx = Mtl(dx,x,ftl,fnl,schtl,schnl,dt,n); 
mxdx = Mnl(xdx,fnl,schnl,dt,n);
mx   = Mnl(x,fnl,schnl,dt,n);

semilogy(abs(alphadx))
figure
semilogy(abs(alphadx2))
figure
semilogy(abs(alphadx3))