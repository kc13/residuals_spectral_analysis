function [S1,f] = runWelch(delta_mat_T,nfft,noverlap,FS,wintype)
%runWelch.m 6/6/22
% note delta_mat assumes time down *rows*, windows in columns

if ~exist('wintype','var')
    wintype = 'hamm';
end

switch wintype
    case 'hamm'
        winarg = nfft;
    case 'rect'
        winarg = ones(nfft,1);
end

data1 = delta_mat_T;
[win_len,n_wins] = size(data1);
zpad = max((nfft-win_len)/2,0);
assert(mod(zpad,2) == 0,'assuming # timepoints / win - nfft would be even'); 

% zero pad windows if short
data1pad = padarray(data1,[zpad,0],0,'both');
tic;
[pow1mat,freq1mat] = arrayfun(@(x) ...
    pwelch(data1pad(:,x),winarg,noverlap,nfft,FS),...
    [1:n_wins],'UniformOutput',false);


S1 = mean(cell2mat(pow1mat),2);
timeElapsed = toc;
f = freq1mat{1};

end %fn

