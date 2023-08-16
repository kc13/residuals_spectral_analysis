function [PDF,lags] = isiPDF_vector(delta_vec)
%isiPDF_vector.m 5/26/22

ISI = getISI_vector(delta_vec);
edges = [1:max(ISI)+1]';
cts = histcounts(ISI,edges)';
PDF = cts./sum(cts);
lags = edges(1:end-1);
end

