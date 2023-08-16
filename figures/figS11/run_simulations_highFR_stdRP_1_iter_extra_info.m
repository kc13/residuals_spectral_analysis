% run_simulations_highFR_stdRP_1_iter_extra_info.m
% to generate a new set of 
% 540 x 1 simulations
% (# param combos x # iter)
% stored out with some extra information


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
};
addpath(scriptdirs{:})
%%
randdir = pwd;
randfile = 'rnd_info_sim_highFR_stdRP_1iter.mat';
randpath = fullfile(randdir,randfile);

write = true; %switch save on/off
outdir = pwd;
%% rnd seed: due to use of parfor, some differences
% in output will remain relative to stored output (see readme)
reuseRndSeed = true;

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
pbase_offset = [35:10:85];  % in Hz
nPB = numel(pbase_offset);
simOpt.run_parfor = false; % parfor for the shuf loops
simOpt.verbose = false;
simOpt.runSplines = false;
nwArr = [30;60;120];
nNW = numel(nwArr);
foArr = [7 9 12 20 32]';
nFO = numel(foArr);
nIter = 1;
fld1arr = {'comp','res'}; % comp = shuffling compensated
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
	    % new suffix so as to not overwrite files already present
        outname = sprintf('len%uf_osc%uiter%u_new.mat',simOpt.len,simOpt.f_osc,nIter);		
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
				mdls = cell(nIter,1); % specific to this script
				rpInfo = cell(nIter,1);
				sim_inputs_all = cell(nIter,1);				
				disp('--------------------------------------------------------------')
				fprintf('running simulations with pbase %g and m %g\n',simOpt.pbase,m);
				disp('--------------------------------------------------------------')
                tic;
                fAll = cell(nIter,1);				
				parfor i = 1:nIter % iteration loop
					fprintf('running iteration %d of %d\n',i,nIter);
					% storing more outputs in this version
					[sim_out,shuf_results,res_results,ISI,sim_inputs] = run_sim_compare_shuf_res(simOpt);
					fAll{i} = toCol(sim_out.f);
					sim_inputs_all{i} = sim_inputs; % specific to this version
					for f1 = 1:nF1
						fld1 = fld1arr{f1};
						S1{i,f1} = toCol(sim_out.(fld1).S1.spec);
						S1sig{i,f1} = toCol(sim_out.(fld1).S1.sigvec);
                        if strcmp(fld1,'res')
                            RP(i) = res_results.res.rpInfo.RPend;
                        end   
						% specific to this version
						if strcmp(fld1,'res')
							mdls{i} = res_results.res.mdls1;
							rpInfo{i} = res_results.res.rpInfo;
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
				sim_data.(fldstrPB).(fldstrM).sim_inputs_all = sim_inputs_all;
				sim_data.(fldstrPB).(fldstrM).res.mdls = mdls;
				sim_data.(fldstrPB).(fldstrM).res.rpInfo = rpInfo;				
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
