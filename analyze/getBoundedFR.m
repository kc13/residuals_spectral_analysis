function [FR,spkPerTr,nPerTr] = getBoundedFR(Recordings,datadir,inds,rt_good,mt_good,boundSpec)
%getBoundedFR.m 7/19/22

startevt = boundSpec.startevt;
endevt = boundSpec.endevt;
start_offset = boundSpec.startoffset;
end_offset = boundSpec.endoffset;
filterFRs = optCheck(boundSpec,'filterFRs',false,[true false]);
min_spk_win = boundSpec.min_spk_win;
win_len_msec = boundSpec.win_len_msec;
R = Recordings;

assert(mod(start_offset,0.001) == 0 && mod(end_offset,0.001) == 0,...
'pre and pst must = integer number of msec')

S = load_spk_evt_data(R,datadir,inds,rt_good,mt_good); 

start_events = notnan(S.Evts.(startevt));
end_events = notnan(S.Evts.(endevt));
assert(numel(start_events) == numel(end_events),...
	'function assumes same number of lower and upper bound events')
	
spk_t = S.units.ts;	

src_raster = interevent_raster(spk_t,start_events,end_events, start_offset, end_offset);

n_trials_raw = size(start_events,1);
assert(size(src_raster,1) == n_trials_raw,...
'trials should be in rows before calling spk_t2delta')


% note: this was originally written to be general 
% enough to accommodate multiple windows per trial
src_delta_raw = cell(n_trials_raw,1);

for t = 1:n_trials_raw
	src_spks = notnan(src_raster(t,:));
	delta_len_msec = (end_events(t)+end_offset - start_events(t)+start_offset)*1000; 
	src_spks(src_spks > delta_len_msec/1000) = [];
	src_delta_full = run_spk_t2delta(src_spks,delta_len_msec);
	src_delta_raw{t} = [src_delta_full];
end %t

n = cellnumel(src_delta_raw);
assert(all(mod(n,win_len_msec) == 0),...
	'src delta not multiple of window length')

src_delta_rs = cellfun(@(x) reshape(x',win_len_msec,[])',...
	src_delta_raw,'UniformOutput',false);
n_wins_raw = size(cell2mat(src_delta_rs),1);

src_nspk = cellsum(src_delta_rs,2);
src_pass = cellgte(src_nspk, min_spk_win);
src_use = cellfun(@(x,ix) x(logical(ix),:),src_delta_rs,src_pass,'UniformOutput',false);

src_delta = cell2mat(src_use);

src_delta_rs = cellfun(@(x) reshape(x',win_len_msec,[])',...
	src_delta_raw,'UniformOutput',false);
n_wins_raw = size(cell2mat(src_delta_rs),1);

src_nspk = cellsum(src_delta_rs,2);
src_pass = cellgte(src_nspk, min_spk_win);
src_use = cellfun(@(x,ix) x(logical(ix),:),src_delta_rs,src_pass,'UniformOutput',false);

src_delta = cell2mat(src_use);

if filterFRs
   src_delta_orig = src_delta;
   [clean_data, info] = find_homogeneous_FR_rows(src_delta);  
   src_delta = arrayfun(@(x) clean_data(x,:), [1:size(clean_data,1)]','UniformOutput',false);  % expected format
end %if

spkPerTr = cellnnz(src_delta);
nPerTr = cellnumel(src_delta); %for fixed 1 ms all 1000s

FR = sum(spkPerTr)/sum(nPerTr)*1000;


end %fn

