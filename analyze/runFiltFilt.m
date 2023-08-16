function [xs] = runFiltFilt(x,smooth_info)
%runFiltFilt.m
% accepts either filter designer object
% or structs w/ coefficients


if isobject(smooth_info)
    % assumes FIR   
    xs = filtfilt(smooth_info.numerator,1,x);
elseif isstruct(smooth_info)
    xs = filtfilt(smooth_info.b,smooth_info.a,x);
else
    error('expected object or struct input')
end
    
end %fn