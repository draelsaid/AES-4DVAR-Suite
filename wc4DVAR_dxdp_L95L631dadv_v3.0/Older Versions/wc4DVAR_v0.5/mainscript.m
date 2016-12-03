clear all
close all

N=40; dt=0.05; T=200;

ftl = @l95tl;
fnl = @l95;
fadj = @l95adj;
dx = zeros(N,T); x  = zeros(N,T); xdx = zeros(N,T);

% Initial conditions
dx(:,1) = 10^(-5)*randn(N,1); x(:,1) = randn(N,1); 
xdx(:,1) = x(:,1)+dx(:,1); 

for k=1:T

[xn,dxn] = Mtl(dx(:,k),x(:,k),ftl,dt); 
[xdxn]   = Mnl(xdx(:,k),fnl,dt);

x(:,k+1) = xn;
dx(:,k+1) = dxn;
xdx(:,k+1) = xdxn;

end
