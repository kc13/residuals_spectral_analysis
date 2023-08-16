function [PDF,lags] = isiPDF2D(delta_mat)
%isiPDF2D.m
% will not look for serial dependencies across rows
ISI = getISI(delta_mat,true);

ixOK = cellnumel(ISI) >= 2;

ISI1row = cellfun(@(x) x(1:end-1), ISI,'UniformOutput',false);
ISI1 = cell2mat(ISI1row(ixOK));
ISI2row = cellfun(@(x) x(2:end), ISI,'UniformOutput',false);
ISI2 = cell2mat(ISI2row(ixOK));

edges = [1:max(union(ISI1,ISI2))+1]'; 

cts = histcounts2(ISI1,ISI2,edges,edges)';
PDF = cts./sum(cts,'all');
lags = edges(1:end-1);

end %fn