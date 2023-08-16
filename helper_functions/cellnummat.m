function [m] = cellnummat(c)
%cellnummat.m 7/6/22

m = cellfun(@(x) str2double(x),c);

end

