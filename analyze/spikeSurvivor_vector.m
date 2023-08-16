function [S,lags] = spikeSurvivor_vector(delta_vec)
%spikeSurvivor_vector.m

ISI = getISI_vector(delta_vec);
edges = [1:max(ISI)+1]';
cts = histcounts(ISI,edges)';
S = cumsum(cts,'reverse')./sum(cts);
lags = edges(1:end-1);

end %fn
