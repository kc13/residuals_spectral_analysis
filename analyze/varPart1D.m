function [init_centers,init_sigma,init_props] = varPart1D(data,k)
%varPart1D.m  12/18/22
% This function is intended to implement the variance partitioning algorithm described by Su and Dy.  See the end of this file for a publication and a Github repository consulted in the development of this function.

    assert(isvector(data),'this function currently only accepts vectors');

    data = toCol(data);
    nO = numel(data);
    labelMat = nan(nO,k);
    labelMat(:,1) = ones(nO,1);

    % sse requires deep learning toolbox
    getSSE = @(x) norm(x-mean(x));

    for i = 1:k-1
        labels = labelMat(:,i);
        SSEs = arrayfun(@(x) getSSE(data(labels == x)),[1:i]);
        % find cluster with max SSE
        maxC = maxind(SSEs,2);
        ixMaxC = labels == maxC;
        maxCmean = mean(data(ixMaxC));
        newLabels = nan(size(labels));
        newLabels(~ixMaxC) = labels(~ixMaxC);
        ixMaxChigh = ixMaxC & (data > maxCmean);
        ixMaxClow = ixMaxC & (data <= maxCmean);
        newLabels(ixMaxClow) = labels(ixMaxClow);
        newLabels(ixMaxChigh) = i+1;
        labelMat(:,i+1) = newLabels;
    end

    finalLabels = labelMat(:,k);

    w = 0;
    init_centers = arrayfun(@(x) mean(data(finalLabels == x)),[1:k]);
    init_sigma = arrayfun(@(x) var(data(finalLabels == x),w),[1:k]); 
    init_props = arrayfun(@(x) mean(finalLabels == x),[1:k]);

end %fn

%%
% References:

% Publication: 
% Su T, Dy JG. In search of deterministic methods for initializing 
% K-means and Gaussian mixture clustering. 
% Intelligent Data Analysis. 2007;11(4):319-38.
% DOI: 10.3233/IDA-2007-11402

% Github reference: 
% https://github.com/tugrulhkarabulut/K-Means-Clustering/
% specific file consulted: initialization_methods.py