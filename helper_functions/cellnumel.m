function [n] = cellnumel(c)
%cellnumel.m 3/13/21

n = cellfun(@(x) numel(x), c);

end

