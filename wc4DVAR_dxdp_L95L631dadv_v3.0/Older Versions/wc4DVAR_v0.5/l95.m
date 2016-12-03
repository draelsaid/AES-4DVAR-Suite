function y=l95(x)
N = length(x); % Number of variables automatically set by user
F=8;
% Periodic Boundary Conditions
    y(1,1) = -x(N-1)*x(N) + x(N)*x(2) - x(1) + F;
    y(2,1) = -x(N)*x(1) + x(1)*x(3) - x(2) + F;
    y(N,1) = -x(N-2)*x(N-1) + x(N-1)*x(1) - x(N) + F;
    
    for i=3:N-1 
         y(i,1) = -x(i-2)*x(i-1) + x(i-1)*x(i+1) - x(i) + F; 
    end
    
end