%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script that fits and plots mean residual waiting-times, hazard rates, and
% generalized SIS compartmental model to gonnorrhea data
%
% Oct 10, 2023 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

%Gonorrhea 2020 incidence data for 15-24 yr olds (obtained from CDC MMWR reports)
incdata = [2111,4111,4109,3864,3342,3513,3130,3212,2824,3036,3381,3114,1614,2316,2044,2990,2562,...
    2686,2907,2492,2952,2818,2898,2863,2914,906,3232,3428,2645,3403,3656,3710,3746,4478,...
    4668,4334,4151,4767,5259,4125,5312,4604,4944,5134,3448,4781,5232,3642,4727,5006,2973,...
    2715,3121];

%parameters for plots
lw= 2; %line width

%time points
t=0:(numel(incdata)-1);

%% Define functions

%define parameter rho as a function of amplitude, frequency, and mean
rho = @(a,w,MU)((1/3)*((36*MU^2*w^2-18*MU^2*w^2*a+1+a+6*sqrt(3)*w*sqrt(16*MU^4*w^4+32*MU^4*w^4*a+16*MU^4*w^4*a^2+8*MU^2*w^2-20*MU^2*w^2*a-MU^2*w^2*a^2+1+a)*MU)*(1+a)^2)^(1/3)/(MU*(1+a))-(1/3)*(12*MU^2*w^2-1)*(1+a)/(MU*((36*MU^2*w^2-18*MU^2*w^2*a+1+a+6*sqrt(3)*w*sqrt(16*MU^4*w^4+32*MU^4*w^4*a+16*MU^4*w^4*a^2+8*MU^2*w^2-20*MU^2*w^2*a-MU^2*w^2*a^2+1+a)*MU)*(1+a)^2)^(1/3))+1/(3*MU));

%mean residual waiting-time (MRWT), derivative of MRWT, and hazard rate
m = @(t,a,w)(rho(a,w,1.91).^2.*(1+a.*cos(w.*t))+4.*w.^2-2.*a.*w.*rho(a,w,1.91).*sin(w.*t))./((1+a.*cos(w.*t)).*rho(a,w,1.91).*(rho(a,w,1.91).^2+4.*w.^2)); 
dm = @(t,a,w)( 2.*a.*w.^2.*(-rho(a,w,1.91).*cos(w.*t)+2.*sin(w.*t).*w-a.*rho(a,w,1.91))./(rho(a,w,1.91).*(rho(a,w,1.91).^2+4.*w.^2).*(1+2.*a.*cos(w.*t)+a.^2.*cos(w.*t).^2)));
eta = @(t,a,w)(-(6.*a.*w.^2.*rho(a,w,1.91).*cos(w.*t)+4.*a.*w.^3.*sin(w.*t)-2.*w.^2.*a.^2.*rho(a,w,1.91)+rho(a,w,1.91).^3+2.*rho(a,w,1.91).^3.*a.*cos(w.*t)+rho(a,w,1.91).^3.*a.^2.*cos(w.*t).^2+4.*rho(a,w,1.91).*w.^2+4.*a.^2.*w.^2.*rho(a,w,1.91).*cos(w.*t).^2)./(-rho(a,w,1.91).^2-2.*rho(a,w,1.91).^2.*a.*cos(w.*t)-rho(a,w,1.91).^2.*a.^2.*cos(w.*t).^2-4.*w.^2-4.*w.^2.*a.*cos(w.*t)+2.*a.*w.*rho(a,w,1.91).*sin(w.*t)+2.*a.^2.*w.*rho(a,w,1.91).*sin(w.*t).*cos(w.*t)));

%duration of infection dist, and eqb dist (survival functions)
P = @(t,x,a,w)(exp(-trapz(0:0.01:t,eta(0:0.01:t,a,w)))/exp(-trapz(0:0.01:x,eta(0:0.01:x,a,w))));
Pe = @(t,x,a,w)(exp(-trapz(0:0.01:t,m(0:0.01:t,a,w)))/exp(-trapz(0:0.01:x,m(0:0.01:x,a,w))));

%check shape of survival functions
%figure(98)
%plot(1:0.1:20,[1 arrayfun(@(t)P(t,1,0.15,2*pi/58),1.1:0.1:20)],1:0.1:20,[1 arrayfun(@(t)Pe(t,1,0.15,2*pi/58),1.1:0.1:20)])

%infectious period dist, and infectious period eqb dist (survival
%functions)
Q = @(t,x,a,w)(exp(-trapz(0:0.01:t,eta(0:0.01:t,a,w))));
Qe = @(t,x,a,w)(exp(-trapz(0:0.01:t,m(0:0.01:t,a,w))));

%check shape of survival functions
%figure(99)
%plot(0:0.1:20,[1 arrayfun(@(t)Q(t,0,0.15,2*pi/58),0.1:0.1:20)],0:0.1:20,[1 arrayfun(@(t)Qe(t,0,0.15,2*pi/58),0.1:0.1:20)])

%integrating factor (mu) (not mean) and the integral of the integrating factor from 0 to t
mu = @(t,a,w,beta)((1+a.*cos(w.*t)).^2.*exp((beta-rho(a,w,1.91)).*t)./(1+a).^2);
intmu = @(t,a,w,beta)((exp(-(-beta+rho(a,w,1.91)).*t)./(beta-rho(a,w,1.91))+4.*a.*((beta-rho(a,w,1.91)).*exp((beta-rho(a,w,1.91)).*t).*cos(w.*t)./((beta-rho(a,w,1.91)).^2+w.^2)+w.*exp((beta-rho(a,w,1.91)).*t).*sin(w.*t)./((beta-rho(a,w,1.91)).^2+w.^2))+a.^2.*((beta-rho(a,w,1.91)).*exp((beta-rho(a,w,1.91)).*t).*cos(2.*w.*t)./((beta-rho(a,w,1.91)).^2+4.*w.^2)+2.*w.*exp((beta-rho(a,w,1.91)).*t).*sin(2.*w.*t)./((beta-rho(a,w,1.91)).^2+4.*w.^2))+a.^2.*exp(-(-beta+rho(a,w,1.91)).*t)./(beta-rho(a,w,1.91)))./(2+4.*a+2.*a.^2)+(1./4).*((1+a).^2.*beta.^4-4.*rho(a,w,1.91).*(1+a).^2.*beta.^3+(3.*((2.*a+2).*rho(a,w,1.91).^2+(a+5./3).*w.^2)).*(1+a).*beta.^2-(6.*(((2./3).*a+2./3).*rho(a,w,1.91).^2+(a+5./3).*w.^2)).*(1+a).*rho(a,w,1.91).*beta+(1+a).^2.*rho(a,w,1.91).^4+(3.*(a+5./3)).*w.^2.*(1+a).*rho(a,w,1.91).^2+2.*w.^4.*(a.^2+2))./((beta.^2-2.*beta.*rho(a,w,1.91)+rho(a,w,1.91).^2+w.^2).*(-beta+rho(a,w,1.91)).*((1./4).*beta.^2-(1./2).*beta.*rho(a,w,1.91)+(1./4).*rho(a,w,1.91).^2+w.^2).*(1+a).^2));

%For large t, mu behaves as:
limitingmu = @(t,a,w,beta)((1+a.*cos(w.*t)).^2./(1+a).^2);

%For large t, integral of mu behaves as:
limitingintmu =@(t,a,w,beta)((1./8).*(-a.^2.*(-beta+rho(a,w,1.91)).^2.*(beta.^2-2.*beta.*rho(a,w,1.91)+rho(a,w,1.91).^2+w.^2).*cos(2.*w.*t)+2.*a.^2.*w.*(-beta+rho(a,w,1.91)).*(beta.^2-2.*beta.*rho(a,w,1.91)+rho(a,w,1.91).^2+w.^2).*sin(2.*w.*t)-(4.*((1./4).*beta.^2-(1./2).*beta.*rho(a,w,1.91)+(1./4).*rho(a,w,1.91).^2+w.^2)).*(4.*a.*(-beta+rho(a,w,1.91)).^2.*cos(w.*t)-4.*a.*w.*(-beta+rho(a,w,1.91)).*sin(w.*t)+(beta.^2-2.*beta.*rho(a,w,1.91)+rho(a,w,1.91).^2+w.^2).*(a.^2+2)))./((-beta+rho(a,w,1.91)).*((1./4).*beta.^2-(1./2).*beta.*rho(a,w,1.91)+(1./4).*rho(a,w,1.91).^2+w.^2).*(beta.^2-2.*beta.*rho(a,w,1.91)+rho(a,w,1.91).^2+w.^2).*(1+a).^2));

%soln curve for gSIS 
i = @(t,a,w,beta,i0)( i0.*mu(t,a,w,beta)./(1+beta.*i0.*intmu(t,a,w,beta)));

%For large t, soln curve for gSIS
istar = @(t,a,w,beta,i0)( limitingmu(t,a,w,beta)./(beta.*limitingintmu(t,a,w,beta)));

%% Define Least square error proc of gSIS to data 
L2 =@(p)((p(5)*p(3)*i(p(6),p(1),p(2),p(3),p(4)).*(1-i(0,p(1),p(2),p(3),p(4)))-incdata(1)).^2 + ...
    trapz((0:52)+p(6),(p(5)*p(3)*i(0:52,p(1),p(2),p(3),p(4)).*(1-i((0:52)+p(6),p(1),p(2),p(3),p(4)))-incdata(1:53)).^2));

%extra options for optimization algorithm
myopt=optimset('Algorithm','interior-point','Display','iter','MaxFunEvals',1500);

%upper and lower bounds for parameters
lb = [-0.999     0  0 0   40000000, 0];
ub = [0.999 0.17 1.535 0.001 45000000, 1000];

% initial cond which alg coverged too
freeparam =[0.157315627770824,0.159642738877891,0.508155568586245,0.000177880447991808,42999721.5035739,97.7774263648788];

%call minimzation alg (for gSIS case)
estparam = fmincon(@(p)L2(p), freeparam,[],[],[],[],lb,ub,[],myopt);

%% Define Least square error proc of SIS to data (features exponential duration of infection distribution)
L2exp =@(p)((p(5)*p(3)*i(p(6),0,0,p(3),p(4)).*(1-i(p(6),0,0,p(3),p(4)))-incdata(1)).^2 + ...
    trapz(0:52,(p(5)*p(3)*i((0:52)+p(6),0,0,p(3),p(4)).*(1-i((0:52)+p(6),0,0,p(3),p(4)))-incdata(1:53)).^2));

%initial cond which alg converged too
freeparam = [0,0.0456197431778813,0.533615365579685,0.000112820516653848,42999712.3111569,3.08317983691195];

%call minimzation alg (for SIS case)
expparam = fmincon(@(p)L2exp(p), freeparam,[],[],[],[],lb,ub,[],myopt);

%%
figure(01)
clf;
%new incidence for gSIS and SIS models
lambda1=estparam(5)*estparam(3).*i(t,estparam(1),estparam(2),estparam(3),estparam(4)).*(1-i(t,estparam(1),estparam(2),estparam(3),estparam(4)));
lambda2=expparam(5)*expparam(3).*i(t,0,0,expparam(3),expparam(4)).*(1-i(t,0,0,expparam(3),expparam(4)));

%plot new incidence from gSIS, SIS, and data
plot(t,incdata,':k',...
    t,lambda1, ...
    t,lambda2, ...
    'linewidth',lw)

xlim([0 max(t)]);
ylabel('New infections')
xlabel('Time (weeks)');
box off
set(gca,'linewidth',lw)
hh=legend('Data','gSIS','SIS');
set(hh,'box','off','location','northwest');

%%
figure(02)
clf;
tau = 0:0.25:300; 
plot(tau,estparam(5)*i(tau,estparam(1),estparam(2),estparam(3),estparam(4))/1000, ...
    tau,estparam(5)*i(tau,0,0,expparam(3),expparam(4))/1000, ...
    'linewidth',lw)
xlim([0 max(tau)]);
ylabel('Incidence (1000s)')
xlabel('Time (weeks)');
box off
set(gca,'linewidth',lw)
hh=legend('gSIS','SIS');
set(hh,'box','off','location','northwest');

%%
figure(03)
clf;
tau = 0:0.5:(2*max(t));
plot(tau,m(tau,estparam(1),estparam(2)), ...
    tau,m(tau,0,0), ...
    'linewidth',lw)
xlim([0 max(tau)]);
ylim([1.8 2.265])
ylabel('Avg. duration of infection (weeks)')
xlabel('Time (weeks)');
box off
set(gca,'linewidth',lw)
hh=legend('gSIS','SIS');
set(hh,'box','off','location','northwest');

%%
figure(04)
clf;
tau = 0:0.25:(2*max(t));
plot(tau,eta(tau,estparam(1),estparam(2)), ...
    tau,eta(tau,0,0), ...
    'linewidth',lw)
xlim([0 max(tau)]);
ylim([0.45 0.57]);
ylabel('\eta(t)')
xlabel('Time (weeks)');
box off
set(gca,'linewidth',lw)
hh=legend('gSIS','SIS');
set(hh,'box','off','location','northwest');
%% Estimation of AIC
% M = # data pts
% k1 = # of free param in gSIS model
% k2 = # of free param in SIS model
% SumSQR = sum of square error

%number of data points (incidence)
M = numel(incdata);

%prevalence of infection
y2=i(t,estparam(1),estparam(2),estparam(3),estparam(4));
baseline=i(t,0,0,expparam(3),expparam(4));

%number of parameters to estimate in gSIS model
k1=5; 

%number of parameters to estimate in SIS model
k2=3; 

%new incidence predicted by gSIS model (formula from above)
%lambda1 = estparam(5)*estparam(3)*y2(1:M).*(1-y2(1:M));

%new incidence predicted by SIS model (formula from above)
%lambda2 = expparam(5)*expparam(3)*baseline(1:M).*(1-baseline(1:M));

%AIC for gSIS model
AICgsis=M*(log(2*pi)+1)+M*log(sum((lambda1-incdata(1:(M))).^2)./M)+2*(k1+1);

%AIC for SIS model
AICsis =M*(log(2*pi)+1)+M*log(sum((lambda2-incdata(1:(M))).^2)./M)+2*(k2+1);