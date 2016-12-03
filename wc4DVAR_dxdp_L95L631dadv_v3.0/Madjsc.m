function dw = Madjsc(dv,x,Madj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjoint model for sc4DVAR                                  %
%   [dw]- Increment output (4D matrix)                       %
%     s - number of model errors (s+1 = number of subwindows)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = size(x,2) - 1;
dw(:,s+1) = dv(:,s+1);

for i = s:-1:1
dw(:,i) = dv(:,i) + Madj(dw(:,i+1),x(:,i));
end

end