%store_delta_matrices_nhp_data.m
% repeats the autocorrelation step only
% taking advantage of the option to save out
% the input delta matrices

%%
clear all;  
close all;
clc;
warning('on','all')

%% add fn paths if not added yet
% assumes script run within local directory
scriptdirs = {
'../../analyze';
'../../helper_functions';
};
addpath(scriptdirs{:})

%%
datadir = '../../spk_evt_data';
addpath(datadir)
tablefile = 'units_table.mat';
tablepath = fullfile(datadir,tablefile);

%% saving prefs
write = true;
outdir = pwd;
outfile = 'delta_mats_GPi_VLa.mat';
outpath = fullfile(outdir,outfile);

%% load Recordings table 
% and ix_lists
fprintf('loading %s\n',tablepath)
load(tablepath)
R = Recordings;

%% initialize variables / set up preferences
area_list = {'GPi';'VLa'}; 
nA = numel(area_list);

%% admitting any reaction/movement time
% since only focused on first second of trial
rt_good = [-Inf Inf];
mt_good = [-Inf Inf];

%% this does an initial screen for a minimum
% number of spikes required for a window to 
% be included. It is set to 0 here, but 
% the filterFRs step will also do post-hoc removal
% of windows with minimal spikes
min_spk_win = 0;

%% specifications for cross-correlation function
% (here, just autocorrelation)
% the function was originally designed to implement post-processing:
% jitter correction and smoothing
% the present manuscript does not make use of these outputs

% for the smoothing filter
% reference:
%https://ocw.mit.edu/courses/9-40-introduction-to-neural-computation-spring-2018/resources/mit9_40s18_lec13/
% (lecture slides from Intro to Neural Computation 
% course taught by Michele Fee)
Fs = 1000;
Fnyq = Fs/2;
cutoff= 300;  
Wn= (cutoff/   Fnyq  );
filtorder = 4; 
[b,a]=butter(filtorder, Wn, 'low'); %Butterworthlow-  pass 
smooth_info.a = a;
smooth_info.b = b;

%% for jitter correction of Amarasingham et al. (2012)
% see reference in interevent_xcorr
nJtr = 400;
jtr_int_msec = 100; 

%% additional preferences for auto-correlation
acorr_prefs.lags_msec = 200;
acorr_prefs.store_input = true;
acorr_prefs.min_spk_win = min_spk_win;
acorr_prefs.filterFRs = true;

%% for setting bounds on trial data used
boundSpec.startevt = 'starttime'; 
boundSpec.startoffset = 0; 
boundSpec.endevt = 'starttime';
boundSpec.endoffset = 1;
% filterFRs = true to filter to trials with roughly 
% homogeneous firing rates
boundSpec.filterFRs = true;
boundSpec.min_spk_win = min_spk_win;
boundSpec.win_len_msec = boundSpec.endoffset*1000;

%%
min_sec = 30; 
min_wins = ceil(min_sec*1000/boundSpec.win_len_msec); 

%%
ixdir = datadir;
ixfile = 'screenedIX.mat';
ixpath = fullfile(ixdir,ixfile);
loadIX = false;
if loadIX
    fprintf('loading %s\n',fullfile(ixdir,ixfile))
    load(ixpath)
else
	for a = 1:nA
		Area = area_list{a};
		minFRdata.(Area) = 1; 
		ix_list = ix_lists.(Area);
		nU = numel(ix_list);
		bdFR.(Area) = nan(nU,1);
		bdNTr.(Area) = nan(nU,1);
		for u = 1:nU
			if mod(u,10)==0; fprintf('screening unit %u of %u from %s\n',u,nU,Area);end  
			[bdFR.(Area)(u),spkTr] = getBoundedFR(R,datadir,ix_list{u},rt_good,mt_good,boundSpec);
			bdNTr.(Area)(u) = length(spkTr);
		end %u
		ix_pass.(Area) = (bdFR.(Area) > minFRdata.(Area) & bdNTr.(Area) >= min_wins);
		ix.(Area) = ix_list(ix_pass.(Area));
	end %a
end
%%
saveIX = false;
if saveIX
	mkdir(ixdir)
	save(ixpath,'ix_lists','ix','ix_pass','bdFR','bdNTr','-v7.3')
end %if
%%

nIX = structfun(@(x) length(x),ix);

for a = 1:nA  
    Area = area_list{a};
    indList = ix.(Area);
    nI = numel(indList);
    u_data = repmat(struct,nI,1);
    clear c_results ev ev_clean 	
	for i = 1:nI
        fprintf('processing %s unit %u of %u\n',Area,i,nI);  
        inds = indList{i}; 	
		S = load_spk_evt_data(R,datadir,inds,rt_good,mt_good); 
		assert(length(S.units) == 1, 'unexpected number of units returned');  
		u_data(i).fname = R.Filename(inds);
		u_data(i).Date = R.Date(inds); 
        u_data(i).Session = R.Session(inds); 
        u_data(i).ix_list_row = i; %row in original ix_list
        u_data(i).chan = S.units.chan;
        u_data(i).sort = S.units.sort;
        u_data(i).brain_area = S.units.brain_area; 
        u_data(i).ri = inds;
		src_spk_t = S.units.ts; 
        evts = S.Evts;
		ev(i) = S.Evts;
        start_events = notnan(S.Evts.(boundSpec.startevt));
        end_events = notnan(S.Evts.(boundSpec.endevt));		
        assert(numel(start_events) == numel(end_events),...
	        'function assumes same number of lower and upper bound events')			
		disp('running autocorrelation')
        c_results(i) = interevent_xcorr(src_spk_t,src_spk_t,...
            start_events,end_events,boundSpec,nJtr,jtr_int_msec,smooth_info,acorr_prefs);		
		clean_ix = c_results(i).cleanIX;
        ev_clean(i) = structfun(@(x) x(clean_ix), ev(i), 'UniformOutput',false);		
	end %i
	unit_data.(Area) = u_data;
    corr_results.(Area) = c_results; 
	evts_clean.(Area) = ev_clean;
end %a

for a = 1:nA
    evtLen.(area_list{a}) = cell2mat(arrayfun(@(y) y.cue_onset-y.starttime, evts_clean.(area_list{a}),'UniformOutput',false)');
end %a
allLens = cell2mat(struct2cell(evtLen));

%%
if write
   fprintf('saving %s\n',outpath)
   save(outpath,'-v7.3') 
end
%%
