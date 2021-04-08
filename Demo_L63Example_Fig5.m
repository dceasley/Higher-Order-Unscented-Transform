%%% Replicates the Figure 5 (d)-(g) of our paper 'A Higher Order Unscented
%%% Transform' 

clear;clc;close all;
addpath('HOUT');

%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%

tau = .1;               %%% Discrete time step
ftime = 10;             %%% Number of model time units to forecast
fsteps=ceil(ftime/tau); %%% Number of forecast steps       

b = 0;          %%% Use the purely deterministic L63, no stochastic forcing

muErr = ones(10,1);
CErr = ones(10,1);
SErr = ones(10,1);
KErr = ones(10,1);
muErrJ = ones(10,1);
CErrJ = ones(10,1);
SErrJ = ones(10,1);
KErrJ = ones(10,1);

runs=500;

for ii=1:runs
    


%%% Choose a random initial condition and run the dynamics to eliminate any transient
x0 = rand(3,1);
[x0,~] = L63(x0,100,tau,b);
x0 = x0(end,:)';

Ens = 10000;
xEns = repmat(x0,1,Ens) + randn(3,Ens)/5;

%%% Using the initial condition x0 which is now on the manifold we produce the data set
N1 = 3;
x1 = L63(xEns,tau*N1,tau,b);
x1 = squeeze(x1(end,:,:));



[mu,C,S,K] = EmpiricalMoments(x1,ones(1,Ens)/Ens);

[sigmas,w] = HigherOrderUnscentedEnsemble(mu,C,S,K,1e-4);
[Jsigmas,Jw] = ScaledUnscentedEnsemble(mu,C);


%%% Ensemble forecast
for i = 1:10
    
    N2 = N1+i;
    x2 = L63(x1,tau*(N2-N1),tau,b);
    x2 = squeeze(x2(end,:,:));
    [muTrue,CTrue,STrue,KTrue] = EmpiricalMoments(x2,ones(1,Ens)/Ens);

    %%% HOUT forecast
    fsigmas = L63(sigmas,tau*(N2-N1),tau,b);
    fsigmas = squeeze(fsigmas(end,:,:));
    [muHOUT,CHOUT,SHOUT,KHOUT] = EmpiricalMoments(fsigmas,w);
    muErr(i) = muErr(i) * norm(muHOUT-muTrue)^(1/runs);
    CErr(i) = CErr(i)  * norm(CHOUT(:)-CTrue(:))^(1/runs);
    SErr(i) = SErr(i) * norm(SHOUT(:)-STrue(:))^(1/runs);
    KErr(i) = KErr(i) * norm(KHOUT(:)-KTrue(:))^(1/runs);

    %%% SUT
    fJsigmas = L63(Jsigmas,tau*(N2-N1),tau,b);
    fJsigmas = squeeze(fJsigmas(end,:,:));
    [muSUT,CSUT,SSUT,KSUT] = EmpiricalMoments(fJsigmas,Jw);
    muErrJ(i) = muErrJ(i) * norm(muSUT-muTrue)^(1/runs);
    CErrJ(i) = CErrJ(i) * norm(CSUT(:)-CTrue(:))^(1/runs);
    SErrJ(i) = SErrJ(i) * norm(SSUT(:)-STrue(:))^(1/runs);
    KErrJ(i) = KErrJ(i) * norm(KSUT(:)-KTrue(:))^(1/runs);

end

end

lg=34;
sm=28;

figure(1);semilogy(muErr,'b','LineWidth',3);hold on;semilogy(muErrJ,'r','LineWidth',3);
set(gca,'fontsize',sm)
ylabel('Error in mean','FontSize',lg)
xlabel('Forecast Steps','FontSize',lg)
xlim([1 10])
h = legend('HOUT','SUT');
set(h,'Location','southeast')
set(h,'fontsize',sm);
print('Figs/errmuL63','-depsc');

figure(2);semilogy(CErr,'b','LineWidth',3);hold on;semilogy(CErrJ,'r','LineWidth',3);
set(gca,'fontsize',sm)
ylabel('Error in variance','FontSize',lg)
xlabel('Forecast Steps','FontSize',lg)
xlim([1 10])
h = legend('HOUT','SUT');
set(h,'Location','southeast')
set(h,'fontsize',sm);
print('Figs/errcovL63','-depsc');

figure(3);semilogy(SErr,'b','LineWidth',3);hold on;semilogy(SErrJ,'r','LineWidth',3);
set(gca,'fontsize',sm)
ylabel('Error in skewness','FontSize',lg)
xlabel('Forecast Steps','FontSize',lg)
xlim([1 10])
h = legend('HOUT','SUT');
set(h,'Location','southeast')
set(h,'fontsize',sm);
print('Figs/errskewL63','-depsc');

figure(4);semilogy(KErr,'b','LineWidth',3);hold on;semilogy(KErrJ,'r','LineWidth',3);
set(gca,'fontsize',sm)
ylabel('Error in kurtosis','FontSize',lg)
xlabel('Forecast Steps','FontSize',lg)
xlim([1 10])
h = legend('HOUT','SUT');
set(h,'Location','southeast')
set(h,'fontsize',sm);
print('Figs/errkurtL63','-depsc');

