% Strong-Constraint 4DVAR

function [DJ,J] = sc4dvar(p,pb,y,fH,fHT,fD,fR,Mnl,Madj)

N = size(pb,1); % Number of variables determined by row size of pb
s = size(pb,2); % Number of subwindows determined by col size of pb
p = p(1:N,1);

xo = p(:,1); xb = pb(:,1);

xmxb = xo - xb;
Binv_xmxb = fD(xmxb);          % Binv*(xo-xb)

Jb = xmxb'*Binv_xmxb;          % Jb term

%[~,Mx] = Mnl(xo);               % M(xo), Non-linear trajectory
Mx = Mnlsc(xo,Mnl,s-1);        % M(xo), Non-linear trajectory

h_M_x = fH(Mx);                % H(M(x)) Obs operator applied to non-lin traj
hMxmy = h_M_x - y;             % Innovations
Rinv_hMxmy = fR(hMxmy);        % R Weighted Innovations

Jo = sum(sum(hMxmy.*Rinv_hMxmy));    % Jo term

J = 0.5*(Jb + Jo);                   % Cost function

% Gradient Calc

% Obtain Adjoint model trajectory
MTHTRhMpmy = Madjsc(fHT(Rinv_hMxmy),Mx,Madj);  
sumMTHTRhMpmy = MTHTRhMpmy(:,1);

DJ = Binv_xmxb + sumMTHTRhMpmy;      % Gradient
end
