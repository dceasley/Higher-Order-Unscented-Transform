function outerproduct = outerProd3(v)
% computes v \otimes v \otimes v

    d = size(v,1);
    outerproduct = repmat(v*v',[1 1 d]).*permute(repmat(v,[1 d d]),[3 2 1]);

end

