function [e_arr] = isemptyCell(arr)
%isemptyCell.m

e_arr = cellfun(@(x) isempty(x), arr);

end

