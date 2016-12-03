function [Cs] = SOAR(nx,Ls,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Second-Order Auto-Regressive correlation function (SOAR)
%
%  Written by Adam El-Said
%
%  List of main variables
%     nx: Numer of gridpoints
%     Ls: Length scale
%      a: radius of circle
% dtheta: Circle angle between each gridpoint
% 
%  Output:
%   [Cg]: Gaussian Correlation Matrix
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = zeros(nx,nx);
Cs = zeros(nx,nx);

dtheta = (2*pi)/nx;
a = L/(2*pi);

for i=1:nx
    for j=1:nx
        theta(i,j) = abs(i-j) * dtheta;
    end
end

for i=1:nx
    for j=1:nx
        Cs(i,j) = ( 1 +  abs(2*a*sin(theta(i,j)/2))  / Ls ) ...
            * exp (   -  abs(2*a*sin(theta(i,j)/2))  / Ls );
    end
end

end

