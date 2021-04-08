function [usigns,utildes,norms] = SymmetricRankOneDecomp4(K,tol)
% This function finds the vectors of a totally symmetric 4-tensor

    residualNorm = norm(K(:));
    
    iter = 1;
    utildes = [];
    while (residualNorm > tol)||(iter==1)
        
        [u, lambda] = HOPM4(K);
        usigns(iter) = sign(lambda);
        utildes(:,iter) = usigns(iter)*nthroot(abs(lambda), 4)*u;
        
        K = K - lambda*outerProd4(u);
        
        residualNorm = norm(K(:));
        norms(iter) = residualNorm;
        
        iter = iter+1;
    end


end

