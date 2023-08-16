function [sim_out,shuf_results,res_results,varargout] = run_sim_compare_shuf_res(simOpt)
%run_sim_compare_shuf_res.m
%creates a simulated spike train
%and creates uncorrected/shuffling/residuals output
%"comp" refers to shuffling "compensation"

% get param settings
len = simOpt.len;
pbase = simOpt.pbase;
nr = simOpt.nr;
k = simOpt.k;
f_osc = simOpt.f_osc;
p_osc = simOpt.p_osc;
osc_st = 0;

disp('running simulation')
[src_spk_delta,src_sim_p,src_sim_p_raw] = simulate_spk(len,pbase,nr,k,f_osc,p_osc,osc_st);
src_spk_t = (find(src_spk_delta)/1000);

stevt = 0;
endevt = len/1000;

[shufopt,resopt] = initpsdopt(simOpt);

% get PSDs
tic; 
shuf_results = one_event_shuf_psd(src_spk_t,stevt,endevt,shufopt);
toc;
tic;
res_results = one_event_res_psd(src_spk_t,stevt,endevt,resopt);
toc;

if resopt.runSplines
    resfld = 'sRes';
else
    resfld = 'res';
end

% extract return values
sim_out.f = toCol(shuf_results.f); % same for res
starr = {shuf_results,res_results};

fld1arr = {'comp',resfld};
fld2arr = {'S1Comp','S1'};

for F = 1:numel(fld1arr)
    st = starr{F}; 
    fld1 = fld1arr{F}; 
    fld2 = fld2arr{F};
    sim_out.(fld1).S1.spec = st.(fld1).(fld2);
    sim_out.(fld1).S1.sigvec = st.(fld1).stats.(fld2).sigvec1TL;
end %F

% three expected output args before varargout
% after: 4th (optional) output is ISIs
if nargout >= 4
    varargout{1} = getISI_vector(src_spk_delta);
end

% 5th (optional) output is the inputs to the simulation
if nargout >= 5
    sim_inputs.src_spk_delta = src_spk_delta;
    sim_inputs.src_sim_p = src_sim_p;
    sim_inputs.src_sim_p_raw = src_sim_p_raw;
    sim_inputs.src_spk_t = src_spk_t;
    varargout{2} = sim_inputs;
end

% 6th (optional) output is autocorrelation
if nargout >= 6
    xcorr_prefs.store_input = true;
    xcorr_prefs.lags_msec = 200;
    varargout{3} = one_event_xcorr(src_spk_t,src_spk_t,stevt,endevt,xcorr_prefs);
end %if



% nested fn
function [shufopt,resopt] = initpsdopt(simOpt)

	% shared between shuf and res
	shufopt.win_len_msec = simOpt.win_len_msec;
	shufopt.min_oscil_spk = 0;
	shufopt.store_stats = true;
	shufopt.store_input = false;
	shufopt.fs_range = [0.1 100];
	shufopt.verbose = simOpt.verbose;
	shufopt.run_parfor = simOpt.run_parfor;
    shufopt.wintype = optCheck(simOpt,'wintype','hamm',{'hamm','rect'});
	
    %copy
    resopt = shufopt;
	
    % now add unique parts
    shufopt.n_shuf = 100;
    shufopt.scope = 'global';
    shufopt.store_shuf = false;	
	
	% resopt-specific
    resopt.min_spk_scope = 'total';
    resopt.runSplines = optCheck(simOpt,'runSplines',false,[true false]);
	if resopt.runSplines
		resopt.splinePrefs = simOpt.splinePrefs;
	end
    resopt.res_options = optCheck(simOpt,'res_options',statset('fitglm'),[]);
    resopt.order = optCheck(simOpt,'order','auto',[]);
    resopt.DispersionFlag = optCheck(simOpt,'DispersionFlag',false,[true false]);
    resopt.distname = optCheck(simOpt,'distname','poisson',[]);
	
end % initpsdopt fn


end %fn