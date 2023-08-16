function [smat] = shuffle_isi_delta_windows(delta_mat)
%shuffle_isi_delta_windows.m
%6/12/22
% assumes nWin x nTime orientation

ISI = getISI(delta_mat);
nI = numel(ISI);
sISI = sort(ISI,'descend');
pool = sISI;
[nW,winLen] = size(delta_mat);
timeLeft = ones(nW,1)*(winLen-1); % -1 to reserve room for init spk
wmat = zeros(nW,winLen);

while ~isempty(pool)
   sz = numel(pool);
   sIX = 1; % was rand
   samp = pool(sIX);
   openIX = find(timeLeft >= samp);
    try
   if numel(openIX) == 1
        oIx = openIX; % randsample treats scalars differently
    else
        oIx = randsample(openIX,1);
   end
    catch
        keyboard
    end
   % now place there, left aligned for now
   if nnz(wmat(oIx,:)) == 0
       wmat(oIx,1) = 1;
   end %if    
   spkPos = find(wmat(oIx,:),1,'last') + samp;
   wmat(oIx,spkPos) = 1;
   timeLeft(oIx) = winLen-spkPos;
   pool(sIX) = []; 
end

tmat = nan(size(wmat));
smat = nan(size(tmat));
for i = 1:nW
    % nothing to shuffle if < 2 spikes but shift will still add some
    % randomization so not left-aligned
    tmat(i,:) = shuffle_isi_delta(wmat(i,:));
    if timeLeft(i) == 0
        shift = 0;
    else
        shift = randsample([0:timeLeft(i)],1);
    end
    smat(i,:) = circshift(tmat(i,:),shift);
end %i


end %fn

