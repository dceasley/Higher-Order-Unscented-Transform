%%% Replicates the Figure 3 and 4 of our paper 'A Higher Order Unscented
%%% Transform' 

clear; clc; close all;
addpath('HOUT');

%%% Generate large ensemble with non-Gaussian distribution
rng(0);
N=100000;
z=randn(2,N); %%% start with Gaussian
A=randn(2,2)/10; B=randn(2,2)/10; C=randn(2,2)/10;
x=A*z+B*z.^2.*sign(z); %%% apply nonlinear function to get non-Gaussian


[mu,C,S,K] = EmpiricalMoments(x,ones(1,N)/N);  %%% Compute moments of the input


%%% Generate HOUT and SUT ensembles
[sigmas,w] = HigherOrderUnscentedEnsemble(mu,C,S,K,1e-5);
[Jsigmas,Jw] = ScaledUnscentedEnsemble(mu,C);


    %%% Plot the initial distribution and HOUT and SUT ensembles
    prob = exp(-sum(z.^2)/2); %%% density of x, inherited from z
    figure(1);scatter(x(1,:),x(2,:),10,prob,'filled');hold on;
    plot(sigmas(1,:),sigmas(2,:),'.r','markersize',30);
    plot(Jsigmas(1,:),Jsigmas(2,:),'.g','markersize',30);
    lg=28;
    sm=22;
    set(gca,'fontsize',sm);
    l=legend('Input Distribution','HOUT','SUT');
    set(l,'fontsize',sm);
    set(l,'Location','southwest')
    xlabel('First Coordinate','fontsize',lg);
    ylabel('Second Coordinate','fontsize',lg);
    print('Figs/initialFig2','-depsc');


a=randn(1,2); b=randn(1,2)/3; %%% Coefficients used in random polynomial
d=length(mu);
cs = .05:0.05:2;    %%% vary strength of the nonlinearities
num = length(cs);

for pow=2:5         %%% run multiple degree polynomials
    
    for i=1:num     %%% vary strength of nonlinearity

        c = cs(i);
        f =@(x) a*x + b*c*x.^pow;   %%% nonlinear polynomial
        %f =@(x) (sqrt(1+c*sum(x.^2))).^pow;  %suggestion by Referee #1
        
        y=f(x);                     %%% apply to ensemble

        muy = mean(y);              %%% estimate true output statistics
        Cy = mean((y-muy).^2);
        Sy = mean((y-muy).^3);
        Ky = mean((y-muy).^4);

        % Estimate statistics with HOUT
        [muyest,Cyest,Syest,Kyest] = EmpiricalMoments(f(sigmas),w);
        
        err_muy_HOUT(i) = abs(muy-muyest);
        err_Cy_HOUT(i) = abs(Cy-Cyest);
        err_Sy_HOUT(i) = abs(Sy-Syest);
        err_Ky_HOUT(i) = abs(Ky-Kyest);
        

        % Estimate statistics with SUT
        [Jmuyest,JCyest,JSyest,JKyest] = EmpiricalMoments(f(Jsigmas),Jw);

        err_muy_SUT(i) = abs(muy-Jmuyest);
        err_Cy_SUT(i) = abs(Cy-JCyest);
        err_Sy_SUT(i) = abs(Sy-JSyest);
        err_Ky_SUT(i) = abs(Ky-JKyest);

    end

    lg=34;
    sm=28;

    figure(2);hold off;
    semilogy(cs,err_muy_HOUT,'b','LineWidth',3);
    hold on;
    semilogy(cs,err_muy_SUT,'r--','LineWidth',3)
    set(gca,'fontsize',sm)
    ylabel('Error in mean','FontSize',lg)
    xlabel('c','FontSize',lg)
    xlim([0 2])
    h = legend('HOUT','SUT');
    set(h,'Location','east')
    set(h,'fontsize',sm);
    print(['Figs/errmu' num2str(pow)],'-depsc');

    figure(3);hold off;
    semilogy(cs,err_Cy_HOUT,'b','LineWidth',3);
    hold on;
    semilogy(cs,err_Cy_SUT,'r--','LineWidth',3)
    set(gca,'fontsize',sm)
    ylabel('Error in variance','FontSize',lg)
    xlabel('c','FontSize',lg)
    h = legend('HOUT','SUT');
    set(h,'Location','east')
    set(h,'fontsize',sm);
    xlim([0 2])
    print(['Figs/errC' num2str(pow)],'-depsc');
    
    figure(4);hold off;
    semilogy(cs,err_Sy_HOUT,'b','LineWidth',3);
    hold on;
    semilogy(cs,err_Sy_SUT,'r--','LineWidth',3)
    set(gca,'fontsize',sm)
    ylabel('Error in skewness','FontSize',lg)
    xlabel('c','FontSize',lg)
    xlim([0 2])
    h = legend('HOUT','SUT');
    set(h,'Location','east')
    set(h,'fontsize',sm);
    print(['Figs/errS' num2str(pow)],'-depsc');
    
    figure(5);hold off;
    semilogy(cs,err_Ky_HOUT,'b','LineWidth',3);
    hold on;
    semilogy(cs,err_Ky_SUT,'r--','LineWidth',3)
    set(gca,'fontsize',sm)
    ylabel('Error in kurtosis','FontSize',lg)
    xlabel('c','FontSize',lg)
    xlim([0 2])
    h = legend('HOUT','SUT');
    set(h,'Location','east')
    set(h,'fontsize',sm);
    print(['Figs/errK' num2str(pow)],'-depsc');

end