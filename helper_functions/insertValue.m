function [newA] = insertValue(A,V,ix)
% assumes vector

    tmpA = toCol(A);

    newA = [tmpA(1:ix-1); V; tmpA(ix:end)];

    if ~iscolumn(A)
        newA = newA';
    end

end

