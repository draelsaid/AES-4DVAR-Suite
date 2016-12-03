% Testing the tangent-linear model with the numerical integration scheme
clear all
close all

addpath( pwd,'L')

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

[Mxdx,Mxdxtraj] = Mtl(dx,x,ftl,fnl,schtl,schnl,dt,n);
[mxdx,mxdxtraj] = Mnl(xdx,fnl,schnl,dt,n);
[mx,mxtraj]   = Mnl(x,fnl,schnl,dt,n);

alphadx(i)  = norm(mxdx - mx)/norm(Mxdx);
alphadx2(i) = alphadx(i)-1;
alphadx3(i) = norm(mxdx - mx - Mxdx);
end

semilogy(abs(alphadx));  figure
semilogy(abs(alphadx2)); figure 
semilogy(abs(alphadx3))

% Validity Test see where the tangent linear model and non-linear model
% begin to diverge away from each other

%This takes a variable and looks at the difference between its trajectory
%and the perturbed trajectory. When the difference is large you will see
%the diff(k) grow 
for k=1:n
p(k)=mxdxtraj(4,k) - mxtraj(4,k);
q(k)=Mxdxtraj(4,k);
diff(k) = p(k) - q(k);
end

% for k=1:n
% plot(mxdxtraj(4,(1:k)) - mxtraj(4,(1:k))); hold; plot(Mxdxtraj(4,(1:k)),'-r')
% drawnow
% end

