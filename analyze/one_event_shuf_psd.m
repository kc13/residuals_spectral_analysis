function [psd_out] = one_event_shuf_psd(src_spk_vec,stevt,endevt,psdopt)
% one_event_shuf_psd.m

win_len_msec = psdopt.win_len_msec;
psd_out.win_len_msec = win_len_msec;

assert(mod((endevt-stevt)*1000,win_len_msec) == 0,...
    'this function assumes that mod(duration,win) = 0')

src_spk = toCol(src_spk_vec);

% historically set at 0
% see below for how scope influences application
% of the min spk param if != 0
min_oscil_spk = psdopt.min_oscil_spk;
psd_out.min_oscil_spk = min_oscil_spk;

% hamming or rectangular window
wintype = optCheck(psdopt,'wintype','hamm',{'hamm','rect'});
psd_out.wintype = wintype;

% shuffling scope
scope = optCheck(psdopt,'scope','global',{'global','localFixed'});
psd_out.scope = scope;

n_shuf = optCheck(psdopt,'n_shuf',0,[]);
psd_out.n_shuf = n_shuf;

store_stats = optCheck(psdopt,'store_stats','false',[true false]);
store_shuf = optCheck(psdopt,'store_shuf','false',[true false]);
store_input = optCheck(psdopt,'store_input','false',[true false]);

% advisable to set parfor to false if alerady parallelism 
% in calling script
run_parfor = optCheck(psdopt,'run_parfor','true',[true false]);

if run_parfor
	poolobj = gcp;
end %if

verbose = optCheck(psdopt,'verbose','true',[true false]);

fs_range = psdopt.fs_range; %inclusive range for sig search

len_msec = (endevt-stevt)*1000;

% note currently inclusive lb and exclusive ub
src_raster = one_event_raster(src_spk,stevt,endevt);

src_delta_vec = spk_t2delta_vector(src_raster,len_msec);

nSpkSrcVec = sum(src_delta_vec);

% if global, min_oscil_spk thresh applies to whole vector
assert(~(strcmp(scope,'global') && (nSpkSrcVec < min_oscil_spk)),...
'fewer than min. required spikes for the spike train');

%trials is synonymous with windows here
n_trials_raw = len_msec/win_len_msec;

% nTr x nT mats
src_delta_rs = reshape(src_delta_vec,win_len_msec,n_trials_raw)';
src_nspk_raw = sum(src_delta_rs,2);

% if local shuf, min oscil spike applied by window
if strcmp(scope,'localFixed') && min_oscil_spk > 0
    src_pass = src_nspk_raw > min_oscil_spike;
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
psd_out.Alpha = 0.05; % is corrected for n points
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

% shuffling correction
if n_shuf > 0
    disp('running shuffle')
    shS1 = nan(nF,n_shuf);
    shf = shS1;
    if store_shuf && store_input
        shufData1 = nan(size(data1,1),n_wins,n_shuf);    
    end
    if run_parfor
        parfor s = 1:n_shuf
            if verbose
                fprintf('sh %d of %d\n',s,n_shuf);   
            end
            sh_src_delta = shuffle_isi_delta_mat(src_delta,2,scope);
            sh_data1 = detrend(sh_src_delta','constant');
            if store_shuf && store_input
                shufData1(:,:,s) = sh_data1;
            end
            [shS1(:,s),shf(:,s)] = runWelch(sh_data1,psd_out.nfft,noverlap,FS,wintype);
        end %s
    else
        for s = 1:n_shuf          
            if verbose
                fprintf('sh %d of %d\n',s,n_shuf);   
            end
            sh_src_delta = shuffle_isi_delta_mat(src_delta,2,scope);
            sh_data1 = detrend(sh_src_delta','constant');
            if store_shuf && store_input
                shufData1(:,:,s) = sh_data1;
            end
            [shS1(:,s),shf(:,s)] = runWelch(sh_data1,psd_out.nfft,noverlap,FS,wintype);
        end %s    
    end %if
    psd_out.comp.mnShS1 = mean(shS1,2);
    psd_out.comp.S1Comp = psd_out.S1./psd_out.comp.mnShS1;
    if store_shuf
        psd_out.comp.shS1 = shS1;
        if store_input
            psd_out.comp.shufData1 = shufData1;
        end
    end
    
    if store_stats
        pow = psd_out.comp.S1Comp;
        % spec stats will do the alpha correction
        [psd_out.comp.stats.S1Comp] = get_spec_stats(pow,psd_out.f,srch_bnds,cntl_bnds,psd_out.Alpha);  
    end
end %n_shuf > 0



end %fn