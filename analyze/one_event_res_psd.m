function [psd_out] = one_event_res_psd(src_spk_vec,stevt,endevt,psdopt)
%one_event_res_psd.m

win_len_msec = psdopt.win_len_msec;
psd_out.win_len_msec = win_len_msec;
assert(mod((endevt-stevt)*1000,win_len_msec) == 0,...
    'this function assumes that mod(duration,win) = 0')

src_spk = toCol(src_spk_vec);

min_oscil_spk = psdopt.min_oscil_spk;

% min spk has been set to 0 for all analyses
psd_out.min_oscil_spk = min_oscil_spk;
% this would only be relevant in the event that min spk > 0
% it determines whether a minimum # of spikes is checked for the whole time series or for the individual windows
min_spk_scope = optCheck(psdopt,'min_spk_scope','total',{'total','wins'});
res_options = optCheck(psdopt,'res_options',statset('fitglm'),[]);
resopt.concatWins = optCheck(psdopt,'concatWins',true,[true false]);

assert(~(resopt.concatWins && isequal(min_spk_scope,'wins')),'concatWins = true && min_spk_scope = wins is not an accepted combination')

resopt.runSplines = optCheck(psdopt,'runSplines',false,[true false]);
if resopt.runSplines
   resopt.splinePrefs = psdopt.splinePrefs;  % will expect method and brkpt
   resfld = 'sRes';
else
   resfld = 'res';
end


wintype = optCheck(psdopt,'wintype','hamm',{'hamm','rect'});
psd_out.wintype = wintype;

% one intercept/RP effect learned for whole time series
resopt.oneModel = true; 

store_stats = optCheck(psdopt,'store_stats','false',[true false]);
store_input = optCheck(psdopt,'store_input','false',[true false]);

% auto uses the RP estimation algorithm
% a number can be input instead if desirable to hard-code
% a specific RP duration
order = optCheck(psdopt,'order','auto',[]);

run_parfor = optCheck(psdopt,'run_parfor','true',[true false]);
verbose = optCheck(psdopt,'verbose','true',[true false]);
resopt.run_parfor = run_parfor;

if run_parfor
    poolobj = gcp;
end

resopt.verbose = verbose; 

fs_range = psdopt.fs_range;  %inclusive range for sig search

% already confirmed integer
len_msec = (endevt-stevt)*1000;

% note currently inclusive lb and exclusive ub
src_raster = one_event_raster(src_spk,stevt,endevt);
src_delta_vec = spk_t2delta_vector(src_raster,len_msec);

nSpkSrcVec = sum(src_delta_vec);
% if checking for a min # of spikes across whole vector, that is done here
assert(nSpkSrcVec >= min_oscil_spk || strcmp('min_spk_scope','wins'),...
'fewer than min. required spikes for the spike train');

%trials synonymous with windows here
n_trials_raw = len_msec/win_len_msec;

% nTr x nT mats
src_delta_rs = reshape(src_delta_vec,win_len_msec,n_trials_raw)';
src_nspk_raw = sum(src_delta_rs,2);

% window removal if local min spk check
if strcmp(min_spk_scope,'wins') && min_oscil_spk > 0
    src_pass = src_nspk_raw > min_oscil_spk;
    src_delta = src_delta_rs(src_pass,:);
else
    src_delta = src_delta_rs;
end

psd_out.n_wins = size(src_delta,1);
psd_out.src_nspk_win = sum(src_delta,2); 
psd_out.srcFR = sum(src_delta(:))/numel(src_delta)*1000;

data1 = detrend(src_delta','constant');

if store_input
    psd_out.data1 = data1;
end

% currently hard-coded significance testing prefs
psd_out.Alpha = 0.05;
FS = 1000;
psd_out.nfft = 2^nextpow2(win_len_msec);
srch_bnds = fs_range;
psd_out.srch_bnds = srch_bnds;
cntl_bnds = [250 500]; % hard coded
psd_out.cntl_bnds = cntl_bnds;

% first: original PSD
disp('running Welch: original spike train')
noverlap = 0;
[psd_out.S1,psd_out.f] = runWelch(data1,psd_out.nfft,noverlap,FS,wintype);

% stats on the original, for reference
if store_stats
    pow = psd_out.S1;
    % spec stats will do the alpha correction
    [psd_out.stats.S1] = get_spec_stats(pow,psd_out.f,srch_bnds,cntl_bnds,psd_out.Alpha);  
end

nF = length(psd_out.f);

% RP estimation (if auto)
resopt.order = order;
[resmat1,mdls1,rpInfoS] = getLastSpkResiduals(src_delta,resopt); 
psd_out.res = struct;  
psd_out.(resfld).mdls1 = mdls1; 
psd_out.(resfld).rpInfo = rpInfoS;

% detrend data in resmat and then nan --> zpad
% operates down cols
rdata1 = detrend(resmat1','constant','omitnan');
rdata1(isnan(rdata1)) = 0;
if store_input
    psd_out.rdata1 = rdata1;
end

disp('running Welch: residuals')
[psd_out.(resfld).S1,psd_out.(resfld).f] = runWelch(rdata1,psd_out.nfft,noverlap,FS,wintype);

% stats on the residuals-corrected PSD
if store_stats
    pow = psd_out.(resfld).S1;
    % spec stats will do the alpha correction
    [psd_out.(resfld).stats.S1] = get_spec_stats(pow,psd_out.f,srch_bnds,cntl_bnds,psd_out.Alpha);  
end

end %fn