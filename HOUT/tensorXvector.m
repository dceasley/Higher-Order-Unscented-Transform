function Stimes1v = tensorXvector(S,v)

    repvec = size(S); %%% for a 3-tensor in R^d, [d d d]
    repvec(1) = 1;    %%% now [1 d d]
                      %%% size(v) is [d 1] so after repmat it is [d d d]
    Stimes1v = squeeze(sum(S.*repmat(v,repvec),1));

end

