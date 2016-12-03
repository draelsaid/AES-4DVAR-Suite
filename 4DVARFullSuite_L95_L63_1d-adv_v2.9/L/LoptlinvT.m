function dw = LoptlinvT(dv,x,Madj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L Operator Adjoint Inverse (linearised) for w-c 4DVAR      %
%   [dw]- Increment output (4D matrix)                       %
%     s - number of model errors (s+1 = number of subwindows)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = size(dv,2) - 1;    % Number of subwindows determined by size of dp
dw(:,s+1) = dv(:,s+1); % First element of vector p is x_0

for i = s:-1:1
dw(:,i) = dv(:,i) + Madj(dw(:,i+1),x(:,i));
end

end