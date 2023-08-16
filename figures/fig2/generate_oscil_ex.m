% generate_oscil_ex.m
% to generate the two oscillation examples
% with same simulation parameters
% as were used for figures 2-3
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
%%
addpath(scriptdirs{:})
datadir = pwd;
datafile = 'figs_2_3_shuf_res_illustrate_data.mat';
datapath = fullfile(datadir,datafile);

%% load file with parameter and random seed info
% (also includes output used for manuscript, 
% but the goal of this script is to generate new output)
vars_to_load = {'H','L'};
D = load(datapath,vars_to_load{:});

%% set to true to recreate the same spike trains that were
% used to generate figures 2-3
% outputs may still differ slightly due to parfor usage
% (see readme)
recreateFigSpikes = false;

%% generate new high and low FR data
nV = numel(vars_to_load);
D2 = struct;
for v = 1:nV
    var = vars_to_load{v};
    if recreateFigSpikes
        rng(D.(var).rngInfo)
    end
    [D2.(var).sim_out,D2.(var).shuf_results,D2.(var).res_results,...
    D2.(var).ISI,D2.(var).sim_inputs,D2.(var).corr_results] = run_sim_compare_shuf_res(D.(var).simOpt);
end %v
%%
H = D2.H;
L = D2.L;
% corrected spectra can be found in:
% shuf_results.comp.S1Comp (comp = shuffling "compensated")
% res_results.res.S1
%%
write = false; % change to save
if write
    outpath = datadir;
    outfile = 'oscil_ex.mat';
    outpath = fullfile(datadir,outfile);
    fprintf('saving %s\n',outpath)
    save(outpath,'-v7.3')
end





