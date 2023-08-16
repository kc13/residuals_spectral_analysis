% run_simulations_pbase_vary.m
% to generate a new set of 
% 144 x 100 simulations
% (# param combos x # iter)
% note: random seeds not available for this case

% (# param combos x # iter)
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
write = true;
outdir = pwd;

%% no random seed option (original random 
% seeds not located)

%% parameters / other settings
simOpt.win_len_msec = 1024;
simOpt.nr = 9;
simOpt.k = 0.7;
simOpt.run_parfor = false; % so can run parfor in outer loop
simOpt.verbose = false;
m_arr = [0:0.2:1]; 
nM = numel(m_arr);
nIter = 100;
veclen = (2^nextpow2(simOpt.win_len_msec))/2+1;
fld1arr = {'comp','res'}; % comp = shuffling compensated
nF1 = numel(fld1arr);
pbase_offset = -1*[2:2:8]; 
nPBO = numel(pbase_offset);

%% outer loops for len, f_osc
nwArr = [30;60;120];
nNW = numel(nwArr);
pbase_arr = [15 20]/1000; 
nPB = numel(pbase_arr); nFO = nPBO;

%%
poolobj = gcp;
%%

for NW = 1:nNW
	nw = nwArr(NW);
	simOpt.len = nw*simOpt.win_len_msec;
	for PB = 1:nPB
        sim_data = struct(); 
        sim_out = struct(); 
		simOpt.pbase = pbase_arr(PB); 
		foArr = simOpt.pbase*1000 + pbase_offset; 
        outname = sprintf('len%upb%uiter%u.mat',simOpt.len,simOpt.pbase*1000,nIter);	
		for FO = 1:nFO
			simOpt.f_osc = foArr(FO);
			fldstrFO = sprintf('fo%g',simOpt.f_osc);	
			for M = 1:nM
				m = m_arr(M);	
				fldstrM = sprintf('m%g',m*100); 
				simOpt.p_osc = m*simOpt.pbase;
				sim_data.(fldstrFO).(fldstrM) = struct;
				sim_data.(fldstrFO).(fldstrM).p_osc = simOpt.p_osc;
                S1 = cell(nIter,nF1);  % for parfor compatibility
                S1sig = cell(nIter,nF1);				
				for f1 = 1:nF1
					fld1 = fld1arr{f1};
                    for i = 1:nIter 
                        S1{i,f1} = nan(veclen,1);
                        S1sig{i,f1} = nan(veclen,1);
                    end %i
				end %f1  
                RP = nan(nIter,1); 
				disp('--------------------------------------------------------------')
				fprintf('running simulations with f_osc %g and m %g\n',simOpt.f_osc,m);
				disp('--------------------------------------------------------------')
                tic;
                fAll = cell(nIter,1);
                parfor i = 1:nIter
                    fprintf('running iteration %d of %d\n',i,nIter);
                    [sim_out,shuf_results,res_results] = run_sim_compare_shuf_res(simOpt);
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
                % unpack
                for f1 = 1:nF1
                    fld1 = fld1arr{f1};
                    for i = 1:nIter
                        sim_data.(fldstrFO).(fldstrM).(fld1).S1(:,i) = S1{i,f1};
                        sim_data.(fldstrFO).(fldstrM).(fld1).S1sig(:,i) = S1sig{i,f1};
                    end %i
                end %f1
                sim_data.(fldstrFO).(fldstrM).res.RP = RP;
                sim_data.(fldstrFO).(fldstrM).f = fAll{1};
                timeElapsed = toc;
                fprintf('time for %u iterations = %u\n',nIter,timeElapsed);	
			end %M
		end	%FO
        % save the file for this setting pair
        if write
            writepath = fullfile(outdir,outname);
            save(writepath,'-v7.3');
        end			
	end %PB
end %NW
