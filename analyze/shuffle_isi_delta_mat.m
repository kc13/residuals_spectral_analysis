function [smat] = shuffle_isi_delta_mat(dmat,dim,scope)
% shuffles the delta vec
% zero regions pre-first/post-final spike left as is

assert(ndims(dmat) <=2,'this fn is not ready for inputs of dim > 2')

if dim == 1
    wmat = dmat';
else
    wmat = dmat;
end

nT = size(wmat,1);

switch scope
    case 'localFixed' % local shuffling
        for t = 1:nT
            omat(t,:) = shuffle_isi_delta(wmat(t,:));
        end
    case 'global' % global shuffling
        wvec = reshape(wmat',[],1);
        ovec = shuffle_isi_delta(wvec);
        omat = reshape(ovec,size(wmat'))';
    case 'windowedGlobal' % window global shuffling
        omat = shuffle_isi_delta_windows(wmat);
    otherwise
        error('unrecognized or not yet supported shuffling scope')
end %sw

if dim == 1
    smat = omat';
else
    smat = omat;
end

end %fn