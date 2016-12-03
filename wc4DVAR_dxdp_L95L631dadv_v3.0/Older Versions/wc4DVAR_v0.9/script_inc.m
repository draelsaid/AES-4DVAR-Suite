% Weak Constraint 4DVAR main sheet
% Adam El-Said, Uni of Reading, 2013

clear all
% Close any open figure windows
g=gcf;
close all
addpath(pwd,'L')

%% Graph position settings
ss=get(0,'ScreenSize'); fig_width=0.55*ss(3); fig_height=0.4*ss(4);
pos_1=[0.02*ss(3),0.5*ss(4),fig_width,fig_height]; pos_2=[0.02*ss(3),0.05*ss(4),fig_width,fig_height]; pos_3=[0.6*ss(3),0.45*ss(4),0.4*ss(3),0.45*ss(4)];
pl=0.6*ss(3); pb=0.05*ss(4); pw=0.3*ss(3); ph=0.35*ss(4);
pos_4=[pl,pb,pw,ph]; pos_5=[0.02*pw,0.8*ph,0.95*pw,0.2*ph]; pos_6=[0.02*pw,0.02*ph,0.95*pw,0.75*ph];

N=40; n=50; dt=0.025;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRUTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnlt = @l95; scht = @rk4nl;
xto = Mnl(rand(N,1),fnlt,scht,dt,n);  % Model 'spin up' to obtain a 'good' true initial state

[~,xttraj] = Mnl(xto,fnlt,scht,dt,n); % TRUE trajectory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model function handle and integration scheme setup
fnl = @l95;      schnl = @rk4nl;
ftl = @l95tl;    schtl = @rk4tl;
fadj= @l95adj;  schadj = @rk4adj;

% These are now the propogators that will be used throughout the
% assimilation, 1 timestep is set since this is how the code works
fnlm = @(x) Mnl(x,fnl,schnl,dt,1);
ftlm = @(dx,x) Mtl(dx,x,ftl,fnl,schtl,schnl,dt,1); 
fadjm = @(dy,y) Madj(dy,y,fadj,fnl,schadj,schnl,dt,1);

% Covariance Matrices (assim)
sigmabt = 0.1; Lsb = 0.1;  % sigma and lengthscale for B
sigmaqt = 0.1; Lsq = 0.2;  % sigma and lengthscale for Q

fcB = @(N) soar(N,Lsb,dt);
fcQ = @(N) soar(N,Lsb,dt);

% Background creation (model error backgrounds (\eta^b) are assumed 0 for now)
fDt = @(x) dmtx(x,sigmabt,sigmaqt);    % D MATRIX function handle
[~,pb] = fDt(xttraj); 
pb = pb + xttraj; 
pb(:,2:end)=0;   % Background x^b_0 is truth perturbed with noise

[~,xbtraj] = Mnl(pb(:,1),fnlt,scht,dt,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Observation creation parameters
sigmaot = 0.001;

fRt = @(x) rmtx(x,sigmaot);      % Setting handle for use inside cost function
[~,y] = rmtx(xttraj,sigmaot);    % Observations are from the truth with some noise
y = y + xttraj;

obs_f_t = 0;      % Obs frequency (time) (Set to 0 if you want to observe ALL times)
obs_f_x = 1;      % Obs frequency (space) (Regular spacing)
obsxo = 1;        % obsxo determines observation starting point (space)
obsto = 1;        % obsto determines observation starting point (time)

fx = floor(N/obs_f_x); % Variable sizes determined by input

fH = @(y) Hop(y,obs_f_x,obs_f_t,obsxo,obsto);         % H operator
fHT = @(y) HopT(y,N,n+1,obs_f_x,obs_f_t,obsxo,obsto); % H operator Transposed

y = fH(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Minimisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigmab = 1; sigmao = 1; sigmaq = 1;

fD = @(x) dmtx(x,sigmab,sigmaq);
fR = @(x) rmtx(x,sigmao);

p = pb; x = pb; dp = zeros(size(p));

b = pb - Lopnl(x,fnlm);
d = y - fH(x);

%% Gradient Test
%   cfdx_inc = @(dx,x) wcdx_inc(dx,x,b,d,fH,fHT,fD,fR,ftlm,fadjm);
%   cfdp_inc = @(dp,p) wcdp_inc(dp,p,b,d,fH,fHT,fD,fR,ftlm,fadjm);
%   phi = testgrad_inc(cfdp_inc,b); 

%% Minimisation parameters
xstr='Maximum number of iterations'; ystr='Tolerance'; zstr={'100','1d-4'};
min_str=inputdlg({xstr,ystr},'Minimization control',1,zstr);
max_iterations=str2num(min_str{1}); tolerance=str2num(min_str{2});

%% Linear Conjugate Gradient
s = size(b,2);
p = reshape(p,s*N,1); dp = reshape(dp,s*N,1);

%% Linear Conjugate gradient (modified matlab pcg giving J on each iteration)
bb = -wcdp_inc(zeros(size(dp)),zeros(size(p)),b,d,fH,fHT,fD,fR,ftlm,fadjm);
afun = @(dp) wcdp_inc(dp,p,b,d,fH,fHT,fD,fR,ftlm,fadjm) + bb;

[X_adp,FLAG,relres,it,dJX,xmat] = minimise_pcg(afun,bb,tolerance,max_iterations,[],[],dp);

% %% MinRes minimisation
% bb = -wcdp_inc(zeros(size(dp)),zeros(size(p)),b,d,fH,fHT,fD,fR,ftlm,fadjm);
% afun = @(dp) wcdp_inc(dp,p,b,d,fH,fHT,fD,fR,ftlm,fadjm) + bb;
% 
% [X_a,FLAG,relres,it,dJX,xmat] = minres(afun,bb,tolerance,max_iterations,[],[],dp);

X_adp = reshape(X_adp,N,s); p = reshape(p,N,s);

pdp = p + X_adp;

xatrajdp = Lopnlinv(pdp,fnlm);
