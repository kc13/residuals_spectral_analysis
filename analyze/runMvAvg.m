function [xs] = runMvAvg(x,smooth_w)
%runMvAvg.m
% 5/22/22 

flt = ones(1,smooth_w)/smooth_w; 
if mod(smooth_w,2) == 0
    xsraw = conv(x,flt);
    xs = xsraw(1:end-1);
else % odd window
    xs = conv(x,flt,'same');
end %if mod


end

