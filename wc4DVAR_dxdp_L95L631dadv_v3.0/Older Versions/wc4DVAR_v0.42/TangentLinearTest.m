clear all
close all

N=40; dt=0.05; T=200;

ftl = @l95tl;
fnl = @l95;
dx = zeros(N,T); x  = zeros(N,T); xdx = zeros(N,T);
alphanorm = zeros(1,16);

for i=1:16
alpha = 10^(1-i);    
dx(:,1) = alpha*randn(N,1); 
x(:,1) = randn(N,1); 
xdx(:,1) = x(:,1) + dx(:,1); 

for k=1:T

dxn = Mtl(dx(:,k),x(:,k),ftl,fnl,dt); 
xdxn   = Mnl(xdx(:,k),fnl,dt);
xn   = Mnl(x(:,k),fnl,dt);

x(:,k+1) = xn;
dx(:,k+1) = dxn;
xdx(:,k+1) = xdxn;

end
alphadx(i)  = norm(xdx(:,T+1) - x(:,T+1))/norm(dx(:,T+1));
alphadx2(i) = alphadx(i)-1;
alphadx3(i) = norm(xdx(:,T+1) - x(:,T+1) - dx(:,T+1));
end
