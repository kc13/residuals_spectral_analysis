function [g] =  cellgte(c,t)
%cellgte.m %3/14/21

g = cellfun(@(x) x >= t, c,'UniformOutput',false);

end

