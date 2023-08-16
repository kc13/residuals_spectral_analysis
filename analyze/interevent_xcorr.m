function [xc_out] = interevent_xcorr(src_spk,targ_spk,start_events,end_events,boundSpec,nJtr,jtr_int_msec,smooth_info,xcorr_prefs)
% interevent_xcorr.m
% this function was originally generated to support
% capabilities beyond basic autocorrelation:
% cross-correlations, jitter correction, smoothing, etc

filterFRs = optCheck(xcorr_prefs,'filterFRs',false,[true false]);
min_spk_win = xcorr_prefs.min_spk_win;
win_len_msec = boundSpec.win_len_msec;

if isfield(xcorr_prefs,'store_input')
    store_input = xcorr_prefs.store_input;
else
    store_input = false;
end

ciLvl = 95;
start_offset = boundSpec.startoffset;
end_offset = boundSpec.endoffset;
lags_msec = xcorr_prefs.lags_msec;

assert(mod(start_offset,0.001) == 0 && mod(end_offset,0.001) == 0,'pre and pst must = integer number of msec');

end_offset_adjusted = end_offset - 0.001;

src_raster = interevent_raster(src_spk, start_events, end_events, start_offset, end_offset);
targ_raster = interevent_raster(targ_spk, start_events, end_events, start_offset, end_offset);

n_trials_raw = size(start_events,1);
assert(size(src_raster,1) == n_trials_raw && size(targ_raster,1) == n_trials_raw,'trials should be in rows before calling spk_t2delta');

src_delta_raw = cell(n_trials_raw,1);
targ_delta_raw = cell(n_trials_raw,1);

for t = 1:n_trials_raw
	src_spks = notnan(src_raster(t,:));
    targ_spks = notnan(targ_raster(t,:));
    delta_len_msec = (end_events(t)+end_offset - start_events(t)+start_offset)*1000; 
    src_spks(src_spks > delta_len_msec/1000) = [];
    targ_spks(targ_spks > delta_len_msec/1000) = [];
    src_delta_full = run_spk_t2delta(src_spks,delta_len_msec);
    targ_delta_full = run_spk_t2delta(targ_spks,delta_len_msec);
    src_delta_raw{t} = src_delta_full;
    targ_delta_raw{t} = targ_delta_full;
end %t

nSDR = cellnumel(src_delta_raw);
assert(all(mod(nSDR,win_len_msec) == 0),...
	'src delta not multiple of window length')
nTDR = cellnumel(targ_delta_raw);
assert(all(mod(nTDR,win_len_msec) == 0),...
	'targ delta not multiple of window length')

src_delta_rs = cellfun(@(x) reshape(x',win_len_msec,[])',...
	src_delta_raw,'UniformOutput',false);
targ_delta_rs = cellfun(@(x) reshape(x',win_len_msec,[])',...
	targ_delta_raw,'UniformOutput',false);
	
n_wins_raw = size(cell2mat(src_delta_rs),1);

src_nspk = cellsum(src_delta_rs,2);
src_pass = cellgte(src_nspk, min_spk_win);

targ_nspk = cellsum(targ_delta_rs,2);
targ_pass = cellgte(targ_nspk, min_spk_win);

isct_pass = celland(src_pass,targ_pass);

src_use = cellfun(@(x,ix) x(logical(ix),:),src_delta_rs,isct_pass,'UniformOutput',false);
targ_use = cellfun(@(x,ix) x(logical(ix),:),targ_delta_rs,isct_pass,'UniformOutput',false);

src_delta_orig = cell2mat(src_use);
targ_delta_orig = cell2mat(targ_use);

if filterFRs
    [src_clean_data,src_info] = find_homogeneous_FR_rows(src_delta_orig);
	[targ_clean_data,targ_info] = find_homogeneous_FR_rows(targ_delta_orig);
    cleanIX = src_info.cleanIX & targ_info.cleanIX;
    src_delta_arr = src_delta_orig(cleanIX,:);
    targ_delta_arr = targ_delta_orig(cleanIX,:);	
    src_delta = arrayfun(@(x) src_delta_arr(x,:), [1:size(src_delta_arr,1)]','UniformOutput',false);
    targ_delta = arrayfun(@(x) targ_delta_arr(x,:), [1:size(targ_delta_arr,1)]','UniformOutput',false);	
    n_trials = size(src_delta,1);
    xc_out.cleanIX = cleanIX;	
end %filterFRs


xc_out.targ_FR = sum([targ_delta{:}])/(numel([targ_delta{:}])/1000);
xc_out.nTargSpk = sum([targ_delta{:}]); % for completeness
xc_out.src_FR = sum([src_delta{:}])/(numel([src_delta{:}])/1000);
xc_out.nSrcSpk = sum([src_delta{:}]);

%pre-allocate
xc_tr = nan(n_trials,(2*lags_msec)+1);
be = nan(1,(2*lags_msec)+1); % will be the same for all
xc_out.n_trials = n_trials;

if store_input
    xc_out.src_delta = src_delta;
    xc_out.targ_delta = targ_delta;
end


[xc,be,xcraw] = getCCF(src_delta,targ_delta,lags_msec);


xc_out.xc = xc;
xc_out.be = be;
xc_out.xcraw = xcraw;

% implements normalization as in Matlab xcorr 'normalized'
xcSraw = sum(cell2mat(arrayfun(@(x) xcorr(src_delta{x,:},src_delta{x,:},lags_msec,'none'),[1:n_trials]','UniformOutput',false)),1);
xcTraw = sum(cell2mat(arrayfun(@(x) xcorr(targ_delta{x,:},targ_delta{x,:},lags_msec,'none'),[1:n_trials]','UniformOutput',false)),1);        
ix0 = lags_msec+1;
normterm = sqrt(xcSraw(ix0)*xcTraw(ix0));
xcn = xcraw./normterm;
xcn2 = despkACF(xcn);

xc_out.xcn = xcn;
xc_out.xcn2 = xcn2;

xc_jtrs = nan(nJtr,size(xc,2));
poolobj = gcp;

% this section implements jitter correction, see ref at end of file
if nJtr > 0
    parfor s = 1:nJtr 
        src_delta_jtr = jitter_delta(src_delta,jtr_int_msec);
        [xc_jtrs(s,:),be] = getCCF(src_delta_jtr,targ_delta,lags_msec);
    end % parfor s
    xc_out.xc_jtrs = xc_jtrs;
    xc_out.xc_jtrs_mn = mean(xc_jtrs,1);
    xc_out.nJtr = nJtr; % as a record
    xc_out.jtr_int_msec = jtr_int_msec;
    xc_corrected = xc-mean(xc_jtrs,1);
    xc_jtrs_corrected = xc_jtrs - mean(xc_jtrs,1);
    xc_out.xc_jtrs_corrected = xc_jtrs_corrected;
    xc_out.xc_corrected = xc_corrected;
    lP = (100-ciLvl)/2;
    uP = lP+ciLvl;
    xc_out.ci_u = prctile(xc_jtrs_corrected,uP,1);
    xc_out.ci_l = prctile(xc_jtrs_corrected,lP,1);
end %if n jtr

% this section implements smoothing, based on either 
% input filter information or moving average method
if isobject(smooth_info)  || isstruct(smooth_info)
	xc_out.smooth_info = smooth_info;
	xc_out.xc_smooth = runFiltFilt(xc,smooth_info);
	if nJtr > 0  % smooth also done after correction
        xc_out.xc_corrected_smooth = runFiltFilt(xc_corrected,smooth_info);		
		xc_jtrs_smooth = nan(size(xc_jtrs));
        for s = 1:nJtr
            xc_jtrs_smooth(s,:) = runFiltFilt(xc_jtrs_corrected(s,:),smooth_info);
        end	
        xc_out.ci_u_smooth = prctile(xc_jtrs_smooth,uP,1);
        xc_out.ci_l_smooth = prctile(xc_jtrs_smooth,lP,1);
		xc_out.xc_jtrs_smooth = xc_jtrs_smooth;
	end %if
elseif smooth_info > 0
	xc_out.smooth_info = smooth_info;
	xc_out.xc_smooth = runMvAvg(xc,smooth_info);
    if nJtr > 0  % smooth done after correction
        xc_out.xc_corrected_smooth = runMvAvg(xc_corrected,smooth_info);
        xc_jtrs_smooth = nan(size(xc_jtrs));
        for s = 1:nJtr
            xc_jtrs_smooth(s,:) = runMvAvg(xc_jtrs_corrected(s,:),smooth_info);
        end
        xc_out.ci_u_smooth = prctile(xc_jtrs_smooth,uP,1);
        xc_out.ci_l_smooth = prctile(xc_jtrs_smooth,lP,1);
		xc_out.xc_jtrs_smooth = xc_jtrs_smooth;
    end		
end %if


end %fn

% jitter reference:
% Amarasingham A, Harrison MT, Hatsopoulos NG, Geman S. Conditional modeling and the jitter method of spike resampling. Journal of Neurophysiology. 2012 Jan 15;107(2):517-31.
