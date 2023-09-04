%build_sim_RP_table_highFR_stdRP.m
% assumes run_simulations_highFR_stdRP.m has been run first
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

%% set write prefs down at bottom

% to construct file names, loops
nIter = 100;
nwArr = [30;60;120];
nNW = numel(nwArr);
simOpt.win_len_msec = 1024;
foArr = [7 9 12 20 32]';
nFO = numel(foArr);
pbase_offset = [35:10:85];
nPB = numel(pbase_offset);
m_arr = [0:0.2:1]; 
nM = numel(m_arr);
fld1arr = {'comp','res'}; 
codes = ['C','R']';
codepairs = nchoosek(codes,2); 
nCO = size(codepairs,1);
fldpairs = nchoosek(fld1arr,2);
nF1 = numel(fld1arr);
% other variables will be
% loaded with the data files

%% table preallocate
tvars = {'NW','FO','PBO','PB','M','I','RP'};
vartypes = {'double','double','double','double','double','double','double'};
nV = numel(tvars);
nR = nNW*nFO*nPB*nM*nIter;
T = table('Size',[nR nV],'VariableNames',tvars,'VariableTypes',vartypes)
iCol = [1:nIter]';

%%
k = 1;
for NW = 1:nNW
    nw = nwArr(NW);
    simOpt.len = nw*simOpt.win_len_msec;
    for FO = 1:nFO
		simOpt.f_osc = foArr(FO);
		f_osc = simOpt.f_osc; 
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
                RP = sim_data.(fldstrPB).(fldstrM).res.RP;
                tblmat = horzcat(repmat([nw f_osc pbo pbSec m],nIter,1),iCol,RP);
                T(k:k+nIter-1,:) = array2table(tblmat);
                k = k+nIter;			
			end %M
		end %pb
	end %FO
end %NW

%%  will write with a suffix to not overwrite existing fig 4 data
write = false;

suffix = '_new';
if write
    outdir = datadir;
    outfile = ['RPtable_highFR_stdRP',suffix,'.mat'];
    outpath = fullfile(outdir,outfile);
    fprintf('saving %s\n',outpath)
    save(outpath,'-v7.3')
end %if
%%
