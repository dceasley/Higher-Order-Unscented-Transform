function [muest,Cest,Sest,Kest] = EmpiricalMoments(sigmas,w)
% calculates the empirical moments

    [d,N] = size(sigmas);

    %Estimated Mean
    muest = sum(repmat(w,d,1).*sigmas,2);

    %Estimated Covariance
    sigmatilde = sigmas-repmat(muest,1,length(sigmas));
    Cest = (sigmatilde.*repmat(w,d,1))*sigmatilde';

    %Estimated Skewness
    Sest = zeros([d d d]);
    for i=1:N
        Sest = Sest + w(i).*outerProd3(sigmas(:,i)-muest);
    end

    %Estimated Kurtosis
    Kest = zeros([d d d d]);
    for i=1:N
        Kest = Kest + w(i).*outerProd4(sigmas(:,i)-muest);
    end

end