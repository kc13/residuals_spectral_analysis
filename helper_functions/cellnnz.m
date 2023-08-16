function [s] = cellnnz(c)
%cellnnz.m %3/13/21

s = cell2mat(cellfun(@(x) nnz(x),c,'UniformOutput',false));

end

