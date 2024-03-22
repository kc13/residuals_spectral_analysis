% fig1_no_oscil_ex.m
% generates figure 1 using stored values used for manuscript
%% 
clear all;
close all;
clc;
%% addpath if not added yet
% assumes script run within local directory
scriptdirs = {'../../plot';'../../helper_functions'};
addpath(scriptdirs{:})
datadir = pwd;
datafile = 'fig1_no_oscil_ex_data.mat';
datapath = fullfile(datadir,datafile);
%%
fprintf('loading %s\n',datapath)
load(datapath)
%%
maxwidin = 7.5; maxhgtin = 8.75;
figure('PaperPosition',[0 0 maxwidin/2 maxhgtin/3.5])
%%
oti = 'synthetic spike train: no oscillation';
otisz = 11;
TL = tlcompact(1,1,oti,otisz);
%%
goArr = gobjects(2,1);
nexttile(TL,1)
figopt = struct();
figopt.doleg = false;
figopt.plotShuf = false;
figopt.xlim = [0 500];
figopt.XT = [0:50:500];
figopt.skipStats = true;
ax = specStatFig(shuf_results,[],'S1',figopt);
tistr = sprintf('base FR = %u Hz',simOpt.pbase*1000);
T = title(tistr);
%% adjust this section to save out .tif
testwrite = true; 
if testwrite
    outdir = pwd;
    outfile = 'fig1_no_oscil_ex.tif';
    outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );
	fprintf('saving %s',strrep(outpath,'.tif','.eps'));
	print(gcf, strrep(outpath,'.tif','.eps'), '-depsc', '-r300' ); 
end
%%
