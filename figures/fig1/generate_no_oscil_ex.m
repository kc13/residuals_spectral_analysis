% generate_no_oscil_ex.m
% to generate no oscillation examples
% with same simulation parameters
% as were used for figure 1
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
datadir = pwd;
datafile = 'fig1_no_oscil_ex_data.mat';
datapath = fullfile(datadir,datafile);
%% load previous parameters
% (not output -- regenerating)
vars_to_load = {'simOpt','rngInfo'};
load(datapath,vars_to_load{:})
%% set to true to recreate the same spike train that was
% used to generate figure 1
% outputs may still differ slightly due to parfor usage
% (see readme)
recreateFig1Spikes = true;
if recreateFig1Spikes
    rng(rngInfo)
end

%%
disp('running simulation and spectral analysis')
[sim_out,shuf_results,res_results,...
    ISI,sim_inputs,corr_results] = run_sim_compare_shuf_res(simOpt);
% uncorrected spectrum can be found in S1 field of both shuf_results and sim_results

%%
write = false; % change to save
if write
    outpath = datadir;
    outfile = 'no_oscil_ex.mat';
    outpath = fullfile(datadir,outfile);
    fprintf('saving %s\n',outpath)
    save(outpath,'-v7.3')
end


