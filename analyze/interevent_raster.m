function [raster] = interevent_raster(spk_t, start_event_t, end_event_t, start_offset, end_offset)
%interevent_raster.m

inbd = @(d,l,u) d(d>=l & d<u);

n_events = length(start_event_t);
assert(n_events == length(end_event_t),...
    'start and end event vector lengths should match');

start_t = start_event_t + start_offset;
end_t = end_event_t + end_offset;

rst = arrayfun(@(st,ed) inbd(spk_t,st,ed)-st, start_t,end_t,'UniformOutput',false);
lens = cellfun(@(x) numel(x), rst);
maxlen = max(lens);

rstmat = nan(n_events,maxlen);
for t = 1:n_events
   if lens(t) > 0
       spks = rst{t};   
       rstmat(t,1:lens(t)) = spks;
   end %if 
end %t

raster = rstmat;

end %fn

