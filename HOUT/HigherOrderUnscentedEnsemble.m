function [sigmas,w] = HigherOrderUnscentedEnsemble(mu,C,S,K,tol,gamma)
% creates sigma points and corresponding weights

    d = length(mu);
    
    vtildes = SymmetricRankOneDecomp3(S,tol/2);
    J = size(vtildes,2);
    
    [signs,utildes] = SymmetricRankOneDecomp4(K,tol/2);
    L = size(utildes,2);
    
    Ctilde = zeros(size(C));
    for i=1:L
        Ctilde = Ctilde + signs(i)*utildes(:,i)*(utildes(:,i))';
    end
    
    lambdaCtildemax = eigs(Ctilde,1); %%% find the largest eigenvalue of Ctilde
    lambdaCmin = eigs(C,1,0);         %%% find smallest eigenvalue of C
    delta = 1.1*sqrt(lambdaCtildemax/lambdaCmin);
    
    while (eigs(C - (1/delta^2)*Ctilde,1,-1)>0)&&(delta>=1)
        delta = delta/2;
    end
    delta = 2*delta;
    
    Chat = C - (1/delta^2)*Ctilde;
    [U,SS,~] = svd(Chat);
    sqrtChat = U*sqrt(SS)*U';
    
    Cbar = zeros(size(K));
    for i=1:d
        Cbar = Cbar + outerProd4(sqrtChat(:,i));
    end
    
    %tol = sqrt(d*norm(Cbar(:)));
    beta = sqrt(tol/2/norm(Cbar(:)));
    
    if (nargin < 6)
        gamma = size(vtildes,2)^(1/3);
    end

    Lhat = sum(signs);
    tildemu = sum(vtildes,2);
    muhat = -gamma^(-2)*tildemu;
    
    muhatO3 = outerProd3(muhat);
    alpha = sqrt(tol/2/norm(muhatO3(:)));
    

    %Sigma Points
    sigmas = zeros(d, 2*d+2*J+2*L+3);

    %sigma_-2 to sigma_0
    sigmas(:,1) = mu;
    sigmas(:,2) = mu + alpha*muhat;
    sigmas(:,3) = mu - alpha*muhat;

    %sigma_1 to sigma_2d
    sigmas(:,4:d+3) = repmat(mu,1,d) + beta*sqrtChat;
    sigmas(:,d+4:2*d+3) = repmat(mu,1,d) - beta*sqrtChat;

    %sigma_{2d+1} to sigma_{2d+2J}
    sigmas(:,2*d+4:2*d+J+3) = repmat(mu,1,J) + gamma*vtildes;
    sigmas(:,2*d+J+4:2*d+2*J+3) = repmat(mu,1,J) - gamma*vtildes;

    %sigma_{2d+2J+1} to sigma_{2d+2J+2L}
    sigmas(:,2*d+2*J+4:2*d+2*J+L+3) = repmat(mu,1,L) + delta*utildes;
    sigmas(:,2*d+2*J+L+4:2*d+2*J+2*L+3) = repmat(mu,1,L) - delta*utildes;


    %Weights
    w = zeros(1,2*d+2*J+2*L+3);
    w(1) = 1-d/beta^2-Lhat/delta^4;
    w(2) = 1/2/alpha;
    w(3) = -1/2/alpha;
    w(4:2*d+3) = 1/(2*beta^2);
    w(2*d+4:2*d+J+3) = 1/(2*gamma^3);
    w(2*d+J+4:2*d+2*J+3) = -1/(2*gamma^3);
    w(2*d+2*J+4:2*d+2*J+L+3) = signs./(2*delta^4);
    w(2*d+2*J+L+4:2*d+2*J+2*L+3) = signs./(2*delta^4);

end