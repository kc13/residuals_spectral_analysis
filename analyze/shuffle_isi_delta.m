function [svec,dshuff,rpIx,isi] = shuffle_isi_delta(dvec)
% shuffles the delta vec
% zero regions pre-first/post-final spike left as is

% 6/6/22 updating to < 2 (either 0 or 1 spk would crash otherwise)
if nnz(dvec) < 2
    svec = dvec;
    warning('< 2 spikes in delta vec, not shuffling')
    return;
end

ixSpk = find(dvec);
dIx = diff(ixSpk);
isi = dIx; 
rpIx = randperm(length(dIx));
dshuff = dIx(rpIx);

svec = zeros(size(dvec));
svec(ixSpk(1)) = 1; 
svec(ixSpk(1)+cumsum(dshuff)) = 1;

% final checks
assert(max(find(svec)) == max(ixSpk),'final ix should match?')
assert(isequal(sort(dIx),sort(diff(find(svec)))),'ISIs should match?')

end %fn