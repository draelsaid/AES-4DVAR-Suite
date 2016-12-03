% Weak Constraint 4DVAR main script
% Adam El-Said, University of Reading (06/2013)
clear all

u=1;uu=1;

% Model settings
N=500; n=2;    % Size of spatial and temporal domains
dt = 0.1;     % 0.025 dt => 3 hour time step (L95), 1 time unit is 5 days
dx = 100/(2*N);  % Spatial interval size
c = -0.3;
C = (c*dt)/dx

fN = N ./(1:ceil(N^0.5)); fN = fN(fN==fix(fN)).' ; fN = unique([fN;N./fN]);
gN = size(fN); obsit = gN(1);

for uu=1:5:21
for u=1:5:21
    
clearvars -except storeitdx storeerrdx storeitdp storeerrdp storeitPdp storeerrPdp ...
    sigmabsigmaq sigmabsigmao sigmaqsigmao storehesssize storekappaSx storekappaPSp ...
    storekappaSp storekappaD storekappaLTDL uu u N dt dx c C fN
n=u;    
% Close any open figure windows
g=gcf;
close all

%% Graph position settings
ss=get(0,'ScreenSize'); fig_width=0.55*ss(3); fig_height=0.4*ss(4);
pos_1=[0.02*ss(3),0.5*ss(4),fig_width,fig_height]; pos_2=[0.02*ss(3),0.05*ss(4),fig_width,fig_height]; pos_3=[0.6*ss(3),0.45*ss(4),0.4*ss(3),0.45*ss(4)];
pl=0.6*ss(3); pbb=0.05*ss(4); pw=0.3*ss(3); ph=0.35*ss(4);
pos_4=[pl,pbb,pw,ph]; pos_5=[0.02*pw,0.8*ph,0.95*pw,0.2*ph]; pos_6=[0.02*pw,0.02*ph,0.95*pw,0.75*ph];

%%%%%%%%%%%%%% REMEMBER TO CHANGE PARAMETERS INSIDE MODEL %%%%%%%%%%%%%%%%%
mlvl=3;

%% TRUTH
fnlt = @adv1d; scht = @lpde;    % TRUE Non-lin model and scheme func handles

% IF THE MODEL IS LINEAR, CHECK MINIMISATION PROCEDURE SUITABILITY.

% True Covariance Matrices D and R
sigmabt = 1; Lsbt = 0.01;  % True sigma_b and lengthscale for C_B
sigmaqt = 1; Lsqt = 0.01;  % True sigma_q and lengthscale for C_Q
sigmaot = 1;               % True sigma_o

fcBt= @(N) soar(N,Lsbt,dx);                    % B MATRIX function handle
fcQt= @(N) laplacian(N,Lsqt,dx);               % Q MATRIX function handle
fDt = @(x) dmtx(x,sigmabt,sigmaqt,fcBt,fcQt);  % D MATRIX function handle
fRt = @(x) rmtx(x,sigmaot);                    % R MATRIX function handle

% Initial Conditions

% Gaussian-shaped curve (for 1D advection equation ONLY)
xc = (N*dx)/2; sigma = 0.1; Ampl = 6; % Mode, Standard dev and Amplitude
xs = linspace(0,1,N);                 % Physical space index 

for i=1:N
    xto(i) = Ampl*exp(-(xs(i)-xc).^2/sigma^2);
end

% % For chaotic ODE models
% xto = Mnl(rand(N,1),fnlt,scht,dt,50); % Model 'spin up' for initial state

% % Truth generation

% Import errors
% etatt=zeros(50,51); etatt(:,2:end)=etat; etat=etatt;

% Generate errors
[~,etat,~,~] =  fDt(rand(N,n+1)); % True model errors
bias = 0; etat = bias + etat;     % Model error bias
etat(:,1) = 0;                    % Model error time adj (eta_0 not exist)

xttraj(:,1) = xto;
for i=1:n
xttraj(:,i+1) = Mnl(xttraj(:,i),fnlt,scht,dt,mlvl) + etat(:,i+1);
end

etat = etat(:,2:end); % Model error time dimension adjustment

% % Time-movie for advection equation
% for i=1:n
%  figure(10), plot(xs,xttraj(:,i),'b'), title('Truth run'), xlabel('Position'), ...
%  legend('Solution'), ylabel('Value')
%  drawnow
% end

%% ASSIMILATION

% Model function handle and integration scheme setup
fnl = @adv1d;      schnl = @lpde;
ftl = @adv1d;      schtl = @lpdetl;
fadj= @adv1d_adj;  schadj = @lpde_adj;

fM = @() adveqn(c,N,dx,dt); % This produces the matrix, soley for condition number calc.

% Weak-constraint model propogators
fnlm = @(x) Mnl(x,fnl,schnl,dt,mlvl);
ftlm = @(dx,x) Mtl(dx,x,ftl,fnl,schtl,schnl,dt,mlvl); 
fadjm = @(dy,y) Madj(dy,y,fadj,fnl,schadj,schnl,dt,mlvl);

% Background error creation [truth plus some noise ~N(0,D)]
[~,xbe,~,~] = fDt(rand(N,n+1));
pb = xttraj + xbe;
pb(:,2:end)=0;             % Background model errors eta^b=0

[~,xbtraj] = Mnl(pb(:,1),fnlt,scht,dt,n);

%% Observation creation parameters

[~,ye,~] = fRt(rand(N,n+1));    
y = xttraj + ye;           % Obs = truth + noise ~N(0,R)

obs_f_t = 1;      % Obs frequency (time) (Regular spacing. =1 obs all times)
obs_f_x = fN(uu);     % Obs frequency (space) (Regular spacing. =1 full obs)
obsxo = 1;        % Determines observation starting point (space)
obsto = 1;        % Determines observation starting point (time)

fx = floor(N/obs_f_x); % Variable sizes determined by input

subwindows = n / obs_f_t

fH = @(y) Hop(y,obs_f_x,obs_f_t,obsxo,obsto);         % H operator
fHT = @(y) HopT(y,N,n+1,obs_f_x,obs_f_t,obsxo,obsto); % H operator Transposed

y = fH(y);

%% Covariance Matrices (ASSIMILATION)

sigmab = 1; Lsb = 0.01;    % Estimated sigmab and lengthscale for B
sigmaq = 1; Lsq = 0.01;   % Estimated sigmaq and lengthscale for Q
sigmao = 1;                % Estimated sigmao

fcB = @(N) soar(N,Lsb,dx);       % C_B MATRIX function handle
fcQ = @(N) laplacian(N,Lsq,dx);  % C_Q MATRIX function handle

fD = @(x) dmtx(x,sigmab,sigmaq,fcB,fcQ);  % D MATRIX function handle
fR = @(x) rmtx(x,sigmao);                       % R MATRIX function handle

%% Gradient Test

% Function handles
       cfdp = @(dp) wcdp(dp,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
       cfdx = @(dx) wcdx(dx,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
  cfsc4dvar = @(dw) sc4dvar(dw,pb,y,fH,fHT,fD,fR,fnlm,fadjm);

% % Tests
%          [phidp,resdp,alpha] = testgrad(cfdp,pb); % Tests wc4dvar dp formulation
%              [phidx,resdx,~] = testgrad(cfdx,pb); % Tests wc4dvar dx formulation
%    [phisc4dvar,ressc4dvar,~] = testgrad(cfsc4dvar,pb,1); % Tests sc4dvar

% % Plot variation of phi with the stepsize alpha
% figure, loglog(alpha,abs(phidp),'b'); hold; loglog(alpha,abs(phidx),'r'); loglog(alpha,abs(phisc4dvar),'g'), ...
% title('Verification of gradient calculation'), xlabel('\alpha'), ylabel('\Phi(\alpha)');

% % Plot variation of residual with the stepsize alpha
% figure, loglog(alpha,resdp,'b'); hold; loglog(alpha,resdx,'r'); loglog(alpha,ressc4dvar,'g'), ...
% title('Variation of residual'), xlabel('\alpha'), ylabel('\Phi(\alpha)-1')

%% Minimisation

% Minimisation parameters
% xstr='Maximum number of iterations'; ystr='Tolerance'; zstr={'1000','1d-4'};
% min_str=inputdlg({xstr,ystr},'Minimization control',1,zstr); 
% max_iterations=str2num(min_str{1}); tolerance=str2num(min_str{2});

max_iterations=2500; tolerance=1d-5;

p = pb; x = pb;            % wc4DVAR Initial guess
xo = p(:,1);               % sc4DVAR initial guess

s = size(pb,2);
p = reshape(p,s*N,1);
x = reshape(x,s*N,1);

%% LINEAR CONJUGATE GRADIENT (requires variable order [DJ,J] in 4dvar)

% % Minimsation for sc4DVAR
% [bbsc,~] = sc4dvar(zeros(size(xo)),pb,y,fH,fHT,fD,fR,fnlm,fadjm); bbsc=-bbsc;
% scfun = @(xo) sc4dvar(xo,pb,y,fH,fHT,fD,fR,fnlm,fadjm) + bbsc;
% 
% [X_asc,FLAGsc,relressc,~,dJXsc,xmatsc] = minimise_pcg(scfun,bbsc,tolerance,max_iterations,[],[],xo);
% 
% % sc4DVAR Solution
% xatrajsc = Mnlsc(X_asc,fnlm,n);

% % Minimisation for Pdp
% Pp = p; 
% fDp = @(Pp) dmtxprecon(Pp,N,n,sigmab,sigmaq,fcB,fcQ); 
% 
% [bbPdp,~] = wcdp(zeros(size(Pp)),pb,y,fH,fHT,fD,fR,fnlm,fadjm); bbPdp = -bbPdp; 
% wcPdpfun = @(Pp) wcdp(Pp,pb,y,fH,fHT,fD,fR,fnlm,fadjm) + bbPdp; 
%  
% [X_aPdp,FLAGPdp,relresPdp,~,dJXPdp,xmatPdp] = ...
%     minimise_pcg(wcPdpfun,bbPdp,tolerance,max_iterations,fDp,[],Pp); 
% 
% % Pdp Solution
% X_aPdp = reshape(X_aPdp,N,s); 
% xatrajPdp = Lopnlinv(X_aPdp,fnlm); 
% Peta = X_aPdp(:,2:end); 
% 
% %%RESET%%
% p = pb; x = pb;            % wc4DVAR Initial guess
% xo = p(:,1);               % sc4DVAR initial guess
% 
% s = size(pb,2);
% p = reshape(p,s*N,1);
% x = reshape(x,s*N,1);
% %%RESET%%

% % Minimsation for dp 
% [bbdp,~] = wcdp(zeros(size(p)),pb,y,fH,fHT,fD,fR,fnlm,fadjm); bbdp=-bbdp; 
% wcdpfun = @(p) wcdp(p,pb,y,fH,fHT,fD,fR,fnlm,fadjm) + bbdp; 
% 
% [X_adp,FLAGdp,relresdp,~,dJXdp,xmatdp] = minimise_pcg(wcdpfun,bbdp,tolerance,max_iterations,[],[],p); 
% 
% % dp Solution 
% X_adp = reshape(X_adp,N,s); 
% xatrajdp = Lopnlinv(X_adp,fnlm); 
% eta = X_adp(:,2:end); 

% Minimisation for dx 
[bbdx,~] = wcdx(zeros(size(x)),pb,y,fH,fHT,fD,fR,fnlm,fadjm); bbdx=-bbdx;
wcdxfun = @(x) wcdx(x,pb,y,fH,fHT,fD,fR,fnlm,fadjm) + bbdx; 

[X_adx,FLAGdx,relresdx,~,dJXdx,xmatdx] = minimise_pcg(wcdxfun,bbdx,tolerance,max_iterations,[],[],x);

% dx Solution
 X_adx = reshape(X_adx,N,s);
 xatrajdx = X_adx;

% Fix iteration count and vector arrangement for plots
% llll = size(dJXPdp); itPdp = llll(1)
% lll = size(dJXdp); itdp = lll(1)
  ll = size(dJXdx); itdx = ll(1)
% l = size(dJXsc); itsc = l(1)
% dJXdp = dJXdp';
% dJXPdp = dJXPdp';
  dJXdx = dJXdx'; 
% dJXsc = dJXsc'; 

% %% MinRes minimisation

% % Minimsation for dp
% [bbdp,~] = wcdp(zeros(size(p)),pb,y,fH,fHT,fD,fR,fnlm,fadjm); bbdp=-bbdp;
%  wcdpfun = @(p) wcdp(p,pb,y,fH,fHT,fD,fR,fnlm,fadjm) + bbdp;
% 
% [X_adp,FLAGdp,relresdp,itdp,dJXdp,xmatdp] = minres(wcdpfun,bbdp,tolerance,max_iterations,[],[],p);

% % dp Solution
% X_adp = reshape(X_adp,N,s);
% xatrajdp = Lopnlinv(X_adp,fnlm);
% eta = X_adp(:,2:end); 

% % Minimsation for dx
% [bbdx,~] = wcdx(zeros(size(x)),pb,y,fH,fHT,fD,fR,fnlm,fadjm); bbdx=-bbdx;
%  wcdxfun = @(x) wcdx(x,pb,y,fH,fHT,fD,fR,fnlm,fadjm) + bbdx;
% 
% [X_adx,FLAGdx,relresdx,itdx,dJXdx,xmatdx] = minres(wcdxfun,bbdx,tolerance,max_iterations,[],[],x);

% % dx Solution
% X_adx = reshape(X_adx,N,s);
% xatrajdx = Lopnlinv(X_adx,fnlm);

%% Compute J (function value) for each iterate (convergence plot)
%  (This is needed for MinRes and min_pcg)
% for i=1:size(xmatdp,2)
%     qq=xmatdp(:,i);
%     [~,Jdp]=wcdp(qq,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
%     JXdp(i)=Jdp;
% end

% for i=1:size(xmatPdp,2)
%     qq=xmatPdp(:,i);
%     [~,JPdp]=wcdp(qq,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
%     JXPdp(i)=JPdp;
% end

for i=1:size(xmatdx,2)
    qq=xmatdx(:,i);
    [~,Jdx]=wcdx(qq,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
    JXdx(i)=Jdx;
end

% for i=1:size(xmatsc,2)
%     qq=xmatsc(:,i);
%     [~,Jsc]=sc4dvar(qq,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
%     JXsc(i)=Jsc;
% end

%% POLAK REBIERRE NON-LINEAR CONJUGATE GRADIENT (requires [J,DJ] in 4dvar
% procedure IN THAT ORDER)

% % Minimisation for dp
% [X_adp JXdp dJXdp itdp] = minimize_mod_crit(p,'wcdp',max_iterations,tolerance,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
% 
% X_adp = reshape(X_adp,N,s);
% xatrajdp = Lopnlinv(X_adp,fnlm);
% 
% % Minimisation for dx
% [X_adx JXdx dJXdx itdx] = minimize_mod_crit(x,'wcdx',max_iterations,tolerance,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
% 
% X_adx = reshape(X_adx,N,s);
% xatrajdx = X_adx;

%% Trajectory time adjustment
% xttraj = xttraj(1:N,2:end);
% xbtraj = xbtraj(1:N,2:end);
% xatrajdx = xatrajdx(1:N,2:end);
% xatrajdp = xatrajdp(1:N,2:end);

%% Error calculations

% Relative total error
terrordx = norm(xttraj - xatrajdx) / norm(xttraj)
%terrordp = norm(xttraj - xatrajdp) / norm(xttraj)
% terrorPdp = norm(xttraj - xatrajPdp) / norm(xttraj)
% terrorsc = norm(xttraj - xatrajsc) / norm(xttraj)

% Relative error of each time interval
% for i=1:s
% errordx(i) = norm(xttraj(:,i) - xatrajdx(:,i))/norm(xttraj(:,i));
% end

% for i=1:s
% errordp(i) = norm(xttraj(:,i) - xatrajdp(:,i))/norm(xttraj(:,i));
% end

% for i=1:s
% errorPdp(i) = norm(xttraj(:,i) - xatrajPdp(:,i))/norm(xttraj(:,i));
% end

% for i=1:s
% errorsc(i) = norm(xttraj(:,i) - xatrajsc(:,i))/norm(xttraj(:,i));
% end

% Relative absolute error spatially and temporally
% for i=1:N
%     for j=1:s
%       cerrordx(i,j) = abs(xttraj(i,j) - xatrajdx(i,j))/abs(xttraj(i,j));
%     end
% end

% for i=1:N
%     for j=1:s
%       cerrorPdp(i,j) = abs(xttraj(i,j) - xatrajPdp(i,j))/abs(xttraj(i,j));
%     end
% end

% for i=1:N
%     for j=1:s
%       cerrordp(i,j) = abs(xttraj(i,j) - xatrajdp(i,j))/abs(xttraj(i,j));
%     end
% end

[kappaSx] = HessNumCondNoCalc(n,N,fD,fR,fM,obs_f_x,mlvl)

storeitdx(u,uu) = itdx;
storeerrdx(u,uu) = terrordx;

%storeitdp(u,uu) = itdp;
%storeerrdp(u,uu) = terrordp;

% storeitPdp(u,uu) = itPdp;
% storeerrPdp(u,uu) = terrorPdp;

 sigmabsigmaq(u,uu) = sigmab/sigmaq;
 sigmabsigmao(u,uu) = sigmab/sigmao;
 sigmaqsigmao(u,uu) = sigmaq/sigmao;
% storehesssize(uu)=(n+1)*N;

storekappaSx(u,uu) = kappaSx;
% storekappaPSp(u,uu) = kappaPSp
% storekappaSp(u,uu) = kappaSp;
%storekappaD(uu) = kappaD;

%storeitdp(storeitdp==0) = []; storeerrdp(storeerrdp==0) = []; storekappaSp(storekappaSp==0) = []; 
%storeitdx(storeitdx==0) = [];  storeerrdx(storeerrdx==0) = []; storekappaSx(storekappaSx==0) = [];  
%sigmabsigmaq(sigmabsigmaq==0) = []; sigmabsigmao(sigmabsigmao==0) = []; sigmaqsigmao(sigmaqsigmao==0) = [];  
%storekappaD(storekappaD==0) = [];

end
end


