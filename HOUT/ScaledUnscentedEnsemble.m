function [sigmas,w] = ScaledUnscentedEnsemble(mu,C,beta)

    d = length(mu);

    if nargin < 3
        beta = sqrt(d);
    end

    sqrtdC = beta*sqrtm(C);
    sigmas = zeros(d,2*d);
    sigmas(:,1:d) = repmat(mu,1,d) + sqrtdC;
    sigmas(:,d+1:2*d) = repmat(mu,1,d) - sqrtdC;
    sigmas(:,2*d+1) = mu;

    w = ones(1,2*d)./(2*beta^2);
    w(2*d+1) = 1-d/beta^2;
    
end