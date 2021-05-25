function [v1, lambda,conv,v2,v3,v4] = HOPM4(K)
% Higher Order Power Method

    % Initial condition based on singular vector of unfolding
    n=length(K);
    Unfolding = reshape(K,n,n^3);
    [u1,~,~]=svds(Unfolding,1);
    if (min(abs(u1))==0)
        u1 = randn(size(u1));
        u1 = u1/norm(u1);
    end
    v1 = u1;
    v2 = u1;
    v3 = u1;
    v4 = u1;

    % HOPM algorithm
    lambda = inf; lamprev=0;
    iter = 1;
    conv = zeros(100,1);
    while (abs(lambda-lamprev) >= 10*eps)&&(iter<100)
        v1 = tensorXvector(tensorXvector(K,v2),v3)*v4;
        conv(iter) = norm(v1);
        v1 = v1/max(eps,norm(v1));
        v2 = tensorXvector(tensorXvector(K,v1),v3)*v4;
        v2 = v2/max(eps,norm(v2));
        v3 = tensorXvector(tensorXvector(K,v1),v2)*v4;
        v3 = v3/max(eps,norm(v3));
        v4 = tensorXvector(tensorXvector(K,v1),v2)*v3;
        v4 = v4/max(eps,norm(v4));
        lamprev=lambda;
        lambda = v1'*(tensorXvector(tensorXvector(K,v1),v1)*v1);
        
        iter = iter+1;
    end
    conv = conv(1:iter-1);

end

