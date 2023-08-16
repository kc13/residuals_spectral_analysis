function [L] = labelPeaks(m,dim)
% labelPeaks.m  10/8/22
% binary output pks or not

assert(dim > 0 && dim < 3,'dim > 2 not yet supported');

    if dim == 1
        mT = m'; % find peaks w/n rows, undo below
    else
        mT = m;
    end
    
    sz = size(mT);
    n = sz(1);

    LT = zeros(sz);
    for r = 1:n
        [~,locs] = findpeaks(mT(r,:));
        LT(r,locs) = 1;
    end %i
    LT = logical(LT);

    if dim == 1
        L = LT'; % return pks down cols
    else
        L = LT;
    end

end