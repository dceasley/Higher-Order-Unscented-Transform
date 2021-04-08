%%% Replicates the Figure 5 (a)-(c) of our paper 'A Higher Order Unscented
%%% Transform' 

clear;clc;close all;
addpath('HOUT');

%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%

tau = .1;               %%% Discrete time step
ftime = 10;             %%% Number of model time units to forecast
fsteps=ceil(ftime/tau); %%% Number of forecast steps       

b = 0;          %%% Use the purely deterministic L63, no stochastic forcing

%%% Choose a random initial condition and run the dynamics to eliminate any transient
seed=3;
rng(seed);
x0 = rand(3,1);
[x0,~] = L63(x0,100,tau,b);
x0 = x0(end,:)';

Ens = 10000;
sigma = 1/5;
xEns = repmat(x0,1,Ens) + sigma*randn(3,Ens);
prob = exp(-sum((xEns-repmat(x0,1,Ens)).^2,1)/(2*sigma^2));


%%% Using the initial condition x0 which is now on the manifold we produce the data set
N1 = 5;
x1 = L63(xEns,tau*N1,tau,b);
x1 = squeeze(x1(end,:,:));


[mu,C,S,K] = EmpiricalMoments(x1,ones(1,Ens)/Ens);

[sigmas,w] = HigherOrderUnscentedEnsemble(mu,C,S,K,1e-2);
[Jsigmas,Jw] = ScaledUnscentedEnsemble(mu,C);


lg=28;
sm=22;


figure(2);
plot(x1(1,:)+x1(2,:),x1(3,:),'.b');hold on;
plot(sigmas(1,:)+sigmas(2,:),sigmas(3,:),'.r','markersize',30);
plot(Jsigmas(1,:)+Jsigmas(2,:),Jsigmas(3,:),'.g','markersize',30);
axis equal;
set(gca,'fontsize',sm);
l=legend('Initial Distribution','HOUT','SUT');
set(l,'fontsize',sm);
set(l,'Location','southwest')
xlabel('x+y','fontsize',lg);
ylabel('z','fontsize',lg);
print(['Figs/initial' num2str(seed)],'-depsc');


figure(1);
x0 = rand(3,1);
[x0,~] = L63(x0,200,tau/10,b);
x0=x0';
x0 = x0(:,5001:end);
plot(x0(1,:)+x0(2,:),x0(3,:),'k');hold on;
plot(x1(1,:)+x1(2,:),x1(3,:),'.b');hold on;
axis equal;



N2=20;
x1 = L63(x1,tau*(N2-N1),tau,b);
x1 = squeeze(x1(end,:,:));
sigmas = L63(sigmas,tau*(N2-N1),tau,b);
sigmas = squeeze(sigmas(end,:,:));
Jsigmas = L63(Jsigmas,tau*(N2-N1),tau,b);
Jsigmas = squeeze(Jsigmas(end,:,:));

figure(1);
plot(x1(1,:)+x1(2,:),x1(3,:),'.r');hold on;
set(gca,'fontsize',sm);
l=legend('Attractor','Initial Distribution','Forecast Distribution');
set(l,'fontsize',sm);
xlabel('x+y','fontsize',lg);
ylabel('z','fontsize',lg);
ylim([0 68]);
print(['Figs/attractor' num2str(seed)],'-depsc');


figure(3);
plot(x1(1,:)+x1(2,:),x1(3,:),'.b');hold on;
plot(sigmas(1,:)+sigmas(2,:),sigmas(3,:),'.r','markersize',30);
plot(Jsigmas(1,:)+Jsigmas(2,:),Jsigmas(3,:),'.g','markersize',30);
axis equal;
set(gca,'fontsize',sm);
l=legend('Forecast Distribution','HOUT','SUT');
set(l,'fontsize',sm);
xlabel('x+y','fontsize',lg);
ylabel('z','fontsize',lg);
print(['Figs/forecast' num2str(seed)],'-depsc');
