function [H,lags] = spikeHazard_vector(delta_vec)
%spikeHazard_vector.m 5/22/22

ISI = getISI_vector(delta_vec);
edges = [1:max(ISI)+1]';
cts = histcounts(ISI,edges)';
PDF = cts./sum(cts);
[S,lags] = spikeSurvivor_vector(delta_vec);
H = PDF./S;


end %fn