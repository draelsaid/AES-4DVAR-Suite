function [M] = adveqn(c,N,dx,dt)

C = c*(dt/dx); % Courant number

r1=zeros(1,N); r1(1)=1+C; r1(2)=-C;
c1=zeros(N,1); c1(1)=1+C; c1(N)=-C;
M=toeplitz(c1,r1);
end
