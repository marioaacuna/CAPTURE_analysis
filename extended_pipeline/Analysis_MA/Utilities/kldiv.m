function d = kldiv(p, q)
    d = sum(p .* log2(p./q));
end