function [PDF,lags] = isiPDF(delta_mat)
%isiPDF.m  5/28/22

ISI = getISI(delta_mat);
edges = [1:max(ISI)+1]';
cts = histcounts(ISI,edges)';
PDF = cts./sum(cts);
lags = edges(1:end-1);

end

