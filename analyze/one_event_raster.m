function [raster] = one_event_raster(spk_t,stevt,endevt)
%one_event_raster.m - 
% bounding and start-subtracting
% 5/22/22

inbd = @(d,l,u) d(d>=l & d<u);

raster = inbd(spk_t,stevt,endevt)-stevt;

end

