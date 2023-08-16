function [s] = cellsum(c,dim)
%cellsum.m %3/13/21

s = cellfun(@(x) sum(x,dim),c,'UniformOutput',false);

end

