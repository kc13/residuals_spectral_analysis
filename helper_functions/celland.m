function [p] = celland(a,b)
%celland.m

p = cellfun(@(x,y) x & y, a,b, 'UniformOutput',false);

end

