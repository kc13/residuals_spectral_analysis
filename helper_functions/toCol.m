function [vc] = toCol(v)
%toCol.m
% vector to col orientation if not already

assert(isvector(v),'this function assumes a vector input')
vc = reshape(v,[],1);

end %fn