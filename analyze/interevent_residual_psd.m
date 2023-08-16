function [psd_out] = interevent_residual_psd(src_spk,start_events,end_events,boundSpec,psdopt)

poolobj = gcp;

resopt.oneModel = optCheck(psdopt,'oneModel',true,[true false]);
resopt.concatWins = false;

win_len_msec = psdopt.win_len_msec;
psd_out.win_len_msec = win_len_msec;

min_oscil_spk = psdopt.min_oscil_spk;
psd_out.min_oscil_spk = min_oscil_spk;
min_spk_win = min_oscil_spk;

resopt.res_options = optCheck(psdopt,'res_options',statset('fitglm'),[]);
% useful for illustrations
storeAllExpMdls = optCheck(psdopt,'storeAllExpMdls',true,[true false]);

% spline tests were only run for synthetic data
% to run on the empirical data, would need to expand the code
% so it can fall back to indicators when the RP duration
% is too short for spline estimation
resopt.runSplines = false;
resfld = 'res';

% order refers to RP duration
order = optCheck(psdopt,'order','auto',[]);

filterFRs = optCheck(psdopt,'filterFRs',false,[true false]);
storeGMMinfo = optCheck(psdopt,'storeGMMinfo',false,[true false]);

store_stats = optCheck(psdopt,'store_stats','false',[true false]);
store_input = optCheck(psdopt,'store_input','false',[true false]);

fs_range = psdopt.fs_range; 

start_offset = boundSpec.startoffset;
end_offset = boundSpec.endoffset;
assert(mod(start_offset,0.001) == 0 && mod(end_offset,0.001) == 0,...
'pre and pst must = integer number of msec');

src_raster = interevent_raster(src_spk, start_events, end_events, start_offset, end_offset);
n_trials_raw = size(start_events,1);
assert(size(src_raster,1) == n_trials_raw,...
'trials should be in rows before calling spk_t2delta');

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

if filterFRs
    src_delta_orig = src_delta;
    [clean_data,info] = find_homogeneous_FR_rows(src_delta);
    src_delta = clean_data;
    if storeGMMinfo
        psd_out.GMMinfo = info;
    end
end

n_wins = size(src_delta,1);
psd_out.n_wins = n_wins;

data1 = detrend(src_delta','constant');

if store_input
    psd_out.data1 = data1;
end

psd_out.Alpha = 0.05;
params.FS = 1000;
psd_out.nfft = 2^nextpow2(win_len_msec);
params.nfft = psd_out.nfft;
srch_bnds = fs_range;
psd_out.srch_bnds = srch_bnds;
cntl_bnds = [250 500]; % hard coded
psd_out.cntl_bnds = cntl_bnds;

disp('running Welch: original spike train')
params.noverlap = 0; 

[psd_out] = run_psd(data1,params,psd_out);

% stats on the uncorrected PSD
if store_stats
	pow = psd_out.S1;
	% spec stats will do the alpha correction
    [psd_out.stats.S1] = get_spec_stats(pow,psd_out.f,srch_bnds,cntl_bnds,psd_out.Alpha);
end %if

nF = length(psd_out.f);


nF = length(psd_out.f);

ISIsrc = getISI(src_delta);
psd_out.ISIsrc = ISIsrc;

resopt.order = psdopt.order;
resopt.run_parfor = true;
resopt.verbose = true;

[resmat1,mdls1,rpInfoS] = getLastSpkResiduals(src_delta,resopt);
psd_out.res = struct;
psd_out.(resfld).mdls1 = mdls1;
if ~storeAllExpMdls
   psd_out.(resfld).rpInfoS = rmfield(rpInfoS,'expMdl'); 
else
   psd_out.(resfld).rpInfoS = rpInfoS;
end

rdata1 = detrend(resmat1','constant','omitnan');
rdata1(isnan(rdata1)) = 0;
if store_input
    psd_out.rdata1 = rdata1;
end

disp('running residuals psd')
[psd_out.(resfld)] = run_psd(rdata1,params,psd_out.(resfld));

% stats on the residuals psd
if store_stats
    pow = psd_out.(resfld).S1;
    % spec stats will do the alpha correction
    [psd_out.(resfld).stats.S1] = get_spec_stats(pow,psd_out.f,srch_bnds,cntl_bnds,psd_out.Alpha);  
end

end %fn

% local fn
function [dataST] = run_psd(data1,params,dataST)
	nfft = params.nfft;
	noverlap = params.noverlap;
    [N,Tr] = size(data1);
	zpad = max((nfft-N)/2,0);
	assert(mod(zpad,2) == 0,...
		'assuming # timepoints / win - nfft would be even'); 
	data1pad = padarray(data1,[zpad,0],0,'both');
	[pow1mat,freq1mat] = arrayfun(@(x) pwelch(data1pad(:,x),...
		nfft,noverlap,nfft,params.FS),[1:Tr],...
		'UniformOutput',false);
	dataST.S1 = mean(cell2mat(pow1mat),2);
	dataST.f = freq1mat{1}; 
end %local fn
