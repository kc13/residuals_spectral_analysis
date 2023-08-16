%build_sim_HR_FA_table_stdFR_RPk0.4nr4.m
% assumes run_simulations_stdFR_RPk0.4nr4.m has been run first
%%
clear all;
close all;
clc;

%% add fn paths if not added yet
% assumes script run within local directory
scriptdirs = {
'../../figures';
};
datadir = pwd;
%%
addpath(scriptdirs{:})
addpath(datadir)
%% set write preferences down towards end

% to construct file names, loops
nIter = 100;
nwArr = [30;60;120];
nNW = numel(nwArr);
simOpt.win_len_msec = 1024;
foArr = [7 9 12 20 32]';
nFO = numel(foArr);
pbase_offset = [2.^(0:5)]; 
nPB = numel(pbase_offset);
m_arr = [0:0.2:1]; 
nM = numel(m_arr);
fld1arr = {'comp','res'}; 
codes = ['C','R']';
codepairs = nchoosek(codes,2); 
nCO = size(codepairs,1);
fldpairs = nchoosek(fld1arr,2);
nF1 = numel(fld1arr);
fvec = [0:500/512:500]';
% other variables will be
% loaded with the data files

% table variable names
A = {'NW','FO','PBO','PB','M'};
B = [join([repmat({'HR'},nF1,1) cellstr(codes)],'')' join([repmat({'FA'},nF1,1) cellstr(codes)],'')'];
C = [join([repmat({'HRD'},nCO,1) cellstr(codepairs)],'')' join([repmat({'FAD'},nCO,1) cellstr(codepairs)],'')'];
tvars = [A,B,C];
nV = numel(tvars);
vartypes = repmat({'double'},nV,1);
nR = nNW*nFO*nPB*nM;
T = table('Size',[nR nV],'VariableNames',tvars,'VariableTypes',vartypes);


%%
k = 1;
for NW = 1:nNW
	nw = nwArr(NW);
	simOpt.len = nw*simOpt.win_len_msec;
	for FO = 1:nFO
		simOpt.f_osc = foArr(FO);
		f_osc = simOpt.f_osc;
		% hit and FA points (m > 0 screen later)
		targIx = knnsearch(fvec,f_osc,'K',3)';
		out5Ix = find(abs(fvec-f_osc)>5);
		pbase_arr = (simOpt.f_osc + pbase_offset)/1000; 
		fname = sprintf('len%uf_osc%uiter%u.mat',simOpt.len,simOpt.f_osc,nIter);
		fpath = fullfile(datadir,fname);
        fprintf('loading %s\n',fpath);
        load(fpath)
		for pb = 1:nPB
			simOpt.pbase = pbase_arr(pb);
			pbase = simOpt.pbase;
			pbSec = simOpt.pbase*1000;
			fldstrPB = sprintf('pb%g',pbSec);
			pbo = pbSec-f_osc;
			for M = 1:nM
				m = m_arr(M);
				fldstrM = sprintf('m%g',m*100); 
				freq = sim_data.(fldstrPB).(fldstrM).f;
				T(k,1:5) = array2table([nw f_osc pbo pbSec m]);
				for f1 = 1:nF1
					fld1 = fld1arr{f1};
					code = codes(f1);
					sig = sim_data.(fldstrPB).(fldstrM).(fld1).S1sig; %nFO x nI
					sigrows = sig(targIx,:);
                    if m ~= 0 
                        hits = ~isnan(sigrows);
                        anyHits = any(hits,1);
                        hitRate = sum(anyHits)/nIter;  
                    else % no hits by definition
                        hitRate = 0;
                    end % if					
					T.(['HR',code])(k) = hitRate;
                    out5rows = sig(out5Ix,:);
                    if m~=0
                        fa = ~isnan(out5rows);
                    else % anything outside 0-100 search range excluded earlier
                        fa = ~isnan(sig); 
                    end %if  
					anyFA = any(fa,1);
					gte1FArate = sum(anyFA)/nIter;%
					T.(['FA',code])(k) = gte1FArate;
				end %f1
				DVs = {'HR','FA'};
				for d = 1:numel(DVs)
					for fp = 1:nCO
						dv = DVs{d};
						fldpair = fldpairs(fp,:);
						fldA = fldpair{1}; fldB = fldpair{2};
						codepair = codepairs(fp,:);
						T.([dv,'D',codepair]) = T.([dv,codepair(1)]) - T.([dv,codepair(2)]);
					end %fp
				end %d
				k = k+1;
			end %M
		end %pb
	end % FO
end %NW

%%
write = false;
suffix = '_new';
if write
    outdir = datadir;
    outfile = ['HRFAtable_stdFR_RPk0.4nr4',suffix,'.mat'];
    outpath = fullfile(outdir,outfile);
    fprintf('saving %s\n',outpath)
    save(outpath,'-v7.3')
end %if
%%
