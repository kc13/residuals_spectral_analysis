function [ISI] = getISI(delta_mat,varargin)
%getISI.m
% 5/7/22


n_trials = size(delta_mat,1);

switch iscell(delta_mat)
    case false
        isiFn = @(x) diff(find(delta_mat(x,:)));
    case true
        isiFn = @(x) diff(find(delta_mat{x,:}));
end

if nargin == 1 || ~varargin{1}
    ISI = cell2mat(arrayfun(@(x) isiFn(x)',[1:n_trials]','UniformOutput',false));
else
    ISI = arrayfun(@(x) isiFn(x)',[1:n_trials]','UniformOutput',false);
end


end %fn