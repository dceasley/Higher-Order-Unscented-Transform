function outerproduct4 = outerProd4(v)
% computes v \otimes v \otimes v \otimes v

    d = size(v,1);
    outerproduct4 = repmat(v,[1 d d d]).*permute(repmat(v,[1 d d d]),[2 1 3 4]).*permute(repmat(v,[1 d d d]),[3 2 1 4]).*permute(repmat(v,[1 d d d]),[2 3 4 1]);

end

