clear all
gamma = .95;
sige = .007
kyratio = 3.2*4;
iyratio = .26;
hbar = .31;
theta = .4;
cyratio = 1 - iyratio;
lbar = 1 - hbar;
A = (1-theta)*lbar/(cyratio*hbar);
delta = iyratio/kyratio;
beta = 1/(theta/kyratio + 1-delta);
kbar = kyratio^(1/(1-theta))*hbar;
ybar = kbar^theta*hbar^(1-theta);
cbar = cyratio*ybar;

aa = [0;0;-kbar;0];
bb = [0;0;(1-delta)*kbar;theta];
cc = [1 -1 1 -1; % coefs on first equation of non-state variables
      hbar lbar 0 0; % coefs on 2nd equation
      0 0 -cbar ybar;
      (1-theta) 0 0 -1];
dd = [0; 0; 0; 1]; % coef on z in each equation

% Stochastic (Euler) Equation
ff = 0;
gg = beta*theta/kyratio;
hh = 0;
jj = [0 0 1 -beta*theta/kyratio];
kk = [0 0 -1 0];
ll = 0;
mm = 0;
nn = gamma;

[pp,qq,rr,ss] = MySOLVE(aa,bb,cc,dd,ff,gg,hh,jj,kk,ll,mm,nn);

% Simulation
nobs = 120;
nsim = 500;
nvars = 6;
hpfilt = hpf(nobs,1600); % converts to H-P filter matrix
randn('state', 159426) % sets seed of random number generator
% don't want to seed according to clock, because we want it to be the same
% each time to make sure changes in model aren't due to different set of
% random variables

% throw away the first 100 observations to lose dependence on initial
% conditions
nm = nobs + 100; 
sim = zeros(nobs,nvars); % initialize matrix

for isim=1:nsim
    r = randn(nm,1); % draw from random normal dist
    r = sige*r;
    z = zeros(nm,1);
    k = zeros(nm,1);
    for iobs=2:nm
        z(iobs,1) = z(iobs-1,1)*nn + r(iobs,1); % populate vector of z's
        k(iobs,1) = pp*k(iobs-1,1) + qq*z(iobs-1,1);
    end
    z = z(101:nm,1);
    k = k(101:nm,1);
    % colon means "all rows" (fill-in first 3 columns)
    sim(:,[1:3]) = (exp(rr([4 1 3],1)*k' + ss([4 1 3],1)*z'))';
    sim(:,5) = exp(k)*kbar;
    sim(:,4) = (exp(pp*k'+qq*z'))'*kbar - (1-delta)*sim(:,5);
    sim(:,6) = sim(:,1)./sim(:,2);
    siml = log(sim);
    [stdpct, corrcon, corrlag] = stat(siml,hpfilt,5);
    stdpctt(:,:,isim) = stdpct;
    corrcont(:,:,isim) = corrcon;
    corrlagt(:,:,isim) = corrlag;
end
disp('Variables are...')
disp('output, hours, consumption, investment, capital, productivity')
statdisp(stdpctt,corrcont,corrlagt)
    
    
    
    
    
    
    
    
    
    
    
    
    