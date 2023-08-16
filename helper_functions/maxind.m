function [I] = maxind(x,dim)
%maxind.m 8/3/21
% for when it is useful to not have to have I be the 2nd output
% inputs: x, dim

[v,I] = max(x,[],dim);

end

