% weak-constraint 4dvar Hessian calculation 
% Adam El-Said uni. reading 21/03/2014

clear all
ti = 40;
%for k=1:ti
    for kk=1:ti
% Cost function parameters
n = 10; s=n;  % Time levels
N = 50;      % Number of gridpoints

% Model parameters 
dx = 0.02;  % Grid spacing
dt = 0.02;  % Time interval
c = -1;    % Advection equation wave speed
  
    tic
    clear D Dinv rtD rtDinv L Linv H

% Cost function matrices
O = zeros(N,N); I = eye(N);

%%%%%%%%%%%%%%%%%%         Correlation matrix D          %%%%%%%%%%%%%%%%%

%%%% LsS=0.1 and LsL=0.16 gives smallest eval in Q and largest eval in B
%%%% for N=500 and dx=0.1
Lsb = 0.002; Lsq = 0.002*kk;    % Correlation lengthscale

% Conditioning ratio settings

sigma_b = 0.1;
sigma_q = 0.05;
sigma_o = 0.05;

  Css = soar(N,Lsb,dx);
%   Cs = soar(N,Lsq,dx);

  [Cl,gamma] = laplacian(N,Lsq,dx);
  
B = sigma_b^2*Css; Binv=inv(B); % B Matrix
%Q = sigma_q^2*Cs;  Qinv=inv(Q); % Q Matrix
Q = sigma_q^2*Cl; Qinv=inv(Q);

% D Matrix
[Dinv{1:s+1,1:s+1}] = deal(O); [D{1:s+1,1:s+1}] = deal(O); [rtDinv{1:s+1,1:s+1}] = deal(O); 
[rtD{1:s+1,1:s+1}] = deal(O); Dinv{1,1} = Binv; D{1,1} = B; rtD{1,1} = B^0.5; 
rtDinv{1,1} = Binv^0.5;

for i=2:s+1
      Dinv{i,i} = Qinv;
      D{i,i} = Q;
      rtD{i,i} = Q^0.5;
      rtDinv{i,i} = Qinv^0.5;
end

Dinv = cell2mat(Dinv);
D = cell2mat(D);
rtDinv = cell2mat(rtDinv);
rtD = cell2mat(rtD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
tic

Mstep = adveqn(c,N,dx,dt); % Model in matrix form (must be LINEAR OR LINEARISED)
M = Mstep^5;

[L{1:s+1,1:s+1}] = deal(O);
for i=1:s+1
    L{i,i} = I;
end

for i=2:s+1
L{i,i-1} = -M;
end

L = cell2mat(L);

[Linv{1:s+1,1:s+1}] = deal(O);
for i=1:s+1
  for j=1:i
    Linv{i,j} = M^(i-j);
  end
end
Linv = cell2mat(Linv);

q = N/2; h=zeros(q,N); j=1;
  for i=1:(N/q):N
    h(j,i) = 1;
    j=j+1;
  end

qt=5; [H{1:n+1,1:n+1}] = deal(zeros(q,N));
  for j=1:qt:n+1
    H{j,j} = h;
  end
H = cell2mat(H);

%% Hessians
G = Linv'*(H'*H)*Linv;
G2 = H*Linv*D*Linv'*H';

Sp = Dinv + sigma_o^(-2)*G;                       % p formulation Hessian
Sx = L'*Dinv*L + sigma_o^(-2)*H'*H;               % x formulation Hessian

%PSp = eye(size(H'*H)) + sigma_o^(-2)*(rtD*G*rtD); % precond p formulation Hessian
%PSpp= eye(size(H*H')) + sigma_o(kk)^(-2)*(G2);        % precond p formulation Hessian (DUAL SPACE)

% Eigenvalue\vector decomposition sorted in order of magnitude
% [VB,DB] = eig(B); DB = sort(diag(DB),'ascend'); % B
% [VQ,DQ] = eig(Q); DQ = sort(diag(DQ),'ascend'); % Q
% [VBi,DBi] = eig(Binv); DBi = sort(diag(DBi),'ascend'); % Binv
% [VQi,DQi] = eig(Qinv); DQi = sort(diag(DQi),'ascend'); % Q
% [VD,DD] = eig(D); DD = sort(diag(DD),'ascend'); % D
% [VG,DG] = eig(G); DG = sort(diag(DG),'ascend'); % Linv^T H^T H Linv
% [VSp,DSp] = eig(Sp); DSp = sort(diag(abs(DSp)),'ascend');
% [VPSp,DPSp] = eig(PSp); DPSp = sort(diag(abs(DPSp)),'ascend');
% [VSx,DSx] = eig(Sx); DSx = sort(diag(abs(DSx)),'ascend');

% [VM,DM] = eig(M); DM = sort(diag(DM),'ascend'); % M

% [VSp,DSp] = eig(Sp); DSp = sort(diag(DSp),'ascend'); % Sp
% [VSx,DSx] = eig(Sx); DSx = sort(diag(DSx),'ascend'); % Sx
% [VpSp,DpSp] = eig(pSp); DpSp = sort(diag(DpSp),'ascend'); % pSp
% [VpSpp,DpSpp] = eig(pSpp); DpSpp = sort(diag(DpSpp),'ascend'); % pSp dual space
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% % Bound calculation
% if abs(DB(1)-DD(1)) < abs(DQ(1)-DD(1))
%     emin = 1;
% else 
%     emin = n;
% end
% 
% if abs(DB(N)-DD(N*(n+1))) < abs(DQ(N)-DD(N*(n+1)))
%     emax = 1;
% else
%     emax = n;
% end

% [psi_min] = psi(DM(1),n,emin); 
% [psi_max] = psi(DM(N),n,emax);
% 
% [gamma_min] = gammma(DM(1),n);
% [omega_min] = omega(DM(1),n);
% [phi_min] = phi(DM(1),n);
% 
% [SumMnorm] = MnormSum(norm(M,inf),n);

% % Unpreconditioned Sp Hessian bounds
% uboundSp(k,kk) = cond(D)*( 1 + [ (DD(1)/sigma_o(kk)^2) * DG(N*(n+1)) ] );
% lboundSp(k,kk) = cond(D)*( ( 1 + (q/(sigma_o(kk)^2*N*emin))*psi_min*DD(1) ) ...
%                   / ( (1 + (q/(sigma_o(kk)^2*N*emax))*psi_max*DD(N*(n+1))) ) );
%     
% % Preconditioned Hessian bounds %% THE SIGMAS ARE NOT IN THE LOWER BOUND
% % BECAUSE lambda_b^2 * DCB(N) = DB(N) %%%
% uboundPSp(k,kk) = 1 + (DD(N*(n+1))/sigma_o(kk)^2)*SumMnorm^2; %1 + (1/sigma_o(kk))*norm(G2,inf);
% lboundPSp(k,kk) = 1 + (q/(N*(n+1)*sigma_o(kk)^2))*[ DB(N)*gamma_min + DQ(N)*omega_min + sqrt(DQ(N))*sqrt(DB(N))*phi_min ];
% 
% % Unpreconditioned Sx Hessian Bounds
% b1max = max(eig(Binv + M'*Qinv*M + sigma_o(kk)^(-2)*h'*h));
% b2max = max(eig(Qinv + M'*Qinv*M + sigma_o(kk)^(-2)*h'*h));
% b3max = max(eig(Qinv + sigma_o(kk)^(-2)*h'*h)); 
% bmaxes = [b1max b2max b3max];
% mmin = min(abs(DM)); mmax = max(abs(DM)); 
% mqmmin = min(eig(M'*Qinv*M)); mqmmax = max(eig(M'*Qinv*M));
% 
% uboundSx(k,kk) = (max(bmaxes) + 2*DQi(end)*real(mmax))/(min(eig(L'*Dinv*L)));
% lboundSx(k,kk) = [(q/(sigma_o(kk)^2*N)) + (n*(DQi(end)+mqmmin) - 2*DQi(end)*real(mmax) + DBi(end))/(n+1)] / ...
%     [(q/(sigma_o(kk)^2*N)) + (n*(DQi(1)+mqmmax) - 2*DQi(1)*real(mmin) + DBi(1))/(n+1)];

% Condition number calculations
kappaD(1,kk) = cond(D);
kappaSp(1,kk) = cond(Sp);
kappaSx(1,kk) = cond(Sx);
%kappaPSp(k,kk) = cond(PSp);
%kappaPSpp(k,kk) = cond(PSpp);
%kappaB(k,kk) = cond(B);
%kappaQ(k,kk) = cond(Q);
%kappaL(k,kk) = cond(L);
%kappaLinv(k,kk) = cond(Linv);
toc

% k loop denotes ROW
% kk loop denotes COLUMN

 % end
% end
    end
%end


% for j=1:ti
%     kappaM(j,1) = cond(M{1,j});
% end