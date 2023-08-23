% run_simulations_stdFR_stdRP_splines.m
% to generate a new set of 
% 540 x 100 simulations
% (# param combos x # iter)
% and create residuals PSDs using
% either indiator functions or splines

%%
clear all;
close all;
clc;
%% add fn paths if not added yet
% assumes script run within local directory
scriptdirs = {
'../../analyze';
'../../figures';
'../../helper_functions';
'../../simulate';
'../../Modified-Spline-Regression';
};
addpath(scriptdirs{:})
%%
randdir = pwd;
randfile = 'rnd_info_sim_stdFR_stdRP_splines.mat';
randpath = fullfile(randdir,randfile);


write = true; %switch save on/off
outdir = pwd;
suffix = '';
%% rnd seed: due to use of parfor, some differences
% in output will remain relative to stored output (see readme)
reuseRndSeed = false;

if reuseRndSeed
	load(randpath)
	rng(rsd,'twister')
end

%% parameters / other settings
simOpt.win_len_msec = 1024;
simOpt.nr = 9;
simOpt.k = 0.7;
m_arr = [0:0.2:1]; 
nM = numel(m_arr);
pbase_offset = [2.^(0:5)]; % in Hz
nPB = numel(pbase_offset);
simOpt.run_parfor = false; % parfor for the shuf loops
simOpt.verbose = false;
splinePrefs.brkpt = 10; 
splinePrefs.s = 0.5;
splinePrefs.method = 'modcard';
splinePrefs.minSplPts = 5;
simOpt.splinePrefs = splinePrefs;

nwArr = [30;60;120];
nNW = numel(nwArr);
foArr = [7 9 12 20 32]';
nFO = numel(foArr);
nIter = 100;
fld1arr = {'res','sRes'};
nF1 = numel(fld1arr);
veclen = (2^nextpow2(simOpt.win_len_msec))/2+1; % length of freq vector

%%
poolobj = gcp; % default is to run parfor over iterations

%%
for NW = 1:nNW
	nw = nwArr(NW);
	simOpt.len = nw*simOpt.win_len_msec;
	for FO = 1:nFO
		sim_data = struct();
		sim_out = struct(); 
		simOpt.f_osc = foArr(FO);
		pbase_arr = (simOpt.f_osc + pbase_offset)/1000; 
        outname = sprintf('len%uf_osc%uiter%u%s.mat',simOpt.len,simOpt.f_osc,nIter,suffix);		
		for pb = 1:nPB
			simOpt.pbase = pbase_arr(pb);
			fldstrPB = sprintf('pb%g',simOpt.pbase*1000);
			for M = 1:nM
				m = m_arr(M);
				simOpt.p_osc = m*simOpt.pbase;
				fldstrM = sprintf('m%g',m*100); 
				sim_data.(fldstrPB).(fldstrM) = struct;
				sim_data.(fldstrPB).(fldstrM).p_osc = simOpt.p_osc;				
				S1 = cell(nIter,nF1);  % temp storage for parfor compatibility
				S1sig = cell(nIter,nF1);
				for f1 = 1:nF1 % preallocation loops
					fld1 = fld1arr{f1};
					for i = 1:nIter 
						S1{i,f1} = nan(veclen,1); 
						S1sig{i,f1} = nan(veclen,1); 
					end %i
				end %f1
				RP = nan(nIter,1); 
				disp('--------------------------------------------------------------')
				fprintf('running simulations with pbase %g and m %g\n',simOpt.pbase,m);
				disp('--------------------------------------------------------------')
                tic;
                fAll = cell(nIter,1);				
				parfor i = 1:nIter % iteration loop
					fprintf('running iteration %d of %d\n',i,nIter);
					[sim_out,res_results,sRes_results] = run_sim_compare_res_sRes(simOpt);
					fAll{i} = toCol(sim_out.f);
					for f1 = 1:nF1
						fld1 = fld1arr{f1};
						S1{i,f1} = toCol(sim_out.(fld1).S1.spec);
						S1sig{i,f1} = toCol(sim_out.(fld1).S1.sigvec);
                        if strcmp(fld1,'res')
                            RP(i) = res_results.res.rpInfo.RPend;
                        end   						
					end %f1
				end %i
				% unpack temp cell arrays
				for f1 = 1:nF1
					fld1 = fld1arr{f1};
					for i = 1:nIter
                        sim_data.(fldstrPB).(fldstrM).(fld1).S1(:,i) = S1{i,f1};
                        sim_data.(fldstrPB).(fldstrM).(fld1).S1sig(:,i) = S1sig{i,f1};					
					end %i
				end %f1
                sim_data.(fldstrPB).(fldstrM).res.RP = RP;
                sim_data.(fldstrPB).(fldstrM).f = fAll{1};
                timeElapsed = toc;				
				fprintf('time for %u iterations = %u\n',nIter,timeElapsed);
			end %M
		end %pb
		% save the file for this len,fosc pair
        if write
            writepath = fullfile(outdir,outname);
            save(writepath,'-v7.3');
        end		
	end %FO
end %NW
