function [vtildes,norms] = SymmetricRankOneDecomp3(S,tol)
% Approximate rank-1 decomposition of a totally symmetric 3-tensor S
    
    residualNorm = norm(S(:));
    
    iter = 1;
    while (residualNorm > tol)||(iter==1)
        
        [v, lambda] = HOPM3(S);
        vtildes(:,iter) = nthroot(lambda, 3)*v;
        
        S = S - lambda*outerProd3(v);
        
        residualNorm = norm(S(:));
        norms(iter) = residualNorm;

        iter = iter+1;
    end

end

