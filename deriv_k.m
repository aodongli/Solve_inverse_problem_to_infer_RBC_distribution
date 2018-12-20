function [res] = deriv_k(dfdk, b, k)
    res = zeros(length(b),1);
    for i = 1:length(b)
        res(i) = dfdk(b(i),k);
    end
end