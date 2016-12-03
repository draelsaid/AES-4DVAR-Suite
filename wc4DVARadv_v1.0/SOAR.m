function [Cs] = SOAR(N,dx,Ls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Second-Order Auto-Regressive correlation function (SOAR)
%
%  Written by Adam El-Said
%
%  List of main variables
%     N: Numer of gridpoints
%     Ls: Length scale
%      a: radius of circle
% dtheta: Circle angle between each gridpoint
% 
%  Output:
%   [Cs]: Second-Order Auto-Regressive Correlation Matrix
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = zeros(N,N);
Cs = zeros(N,N);
dtheta = (2*pi)/N;
a = (N*dx)/(2*pi);

for i=1:N
    for j=1:N
        theta(i,j) = abs(i-j) * dtheta;
    end
end

for i=1:N
    for j=1:N
        Cs(i,j) = ( 1 +  abs(2*a*sin(theta(i,j)/2))  / Ls ) ...
            * exp (   -  abs(2*a*sin(theta(i,j)/2))  / Ls );
    end
end

end

