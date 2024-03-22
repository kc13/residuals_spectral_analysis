% figS1_RP_no_oscil.m
% generates Fig. S1 using stored values

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
colorfile = 'priColors.mat';
colordir = '../../plot';
colorpath = fullfile(colordir,colorfile);

%%
fprintf('loading %s\n',datapath)
N = load(datapath);

%%
fprintf('loading %s\n',colorpath)
load(colorpath);

%% prefs
xLshort = [1 20];
xLmore = [1 200];
xlbl = sprintf('ISI: zoom to [%u, %u] ms',xLmore);
rpInfo = N.res_results.res.rpInfoS;
HLix = [1 5 rpInfo.optLag];

maxwidin = 7.5; maxhgtin = 8.75; 
figure('PaperPosition',[0 0 maxwidin/2 maxhgtin/2])

oti = 'synthetic spike train: no oscillation';
otisz = 11;
TL = tlcompact(2,1,oti,otisz);

%%
figprefs = struct();
simOpt = N.simOpt;
m = N.m;
nexttile(TL) %A
figprefs.xlims = xLmore;
figprefs.xlbl = xlbl;
figprefs.doleg = true;
plotExpToISI(rpInfo,HLix,priColors,figprefs);
tistr = sprintf('base FR = %u Hz',simOpt.pbase*1000);
title(tistr)
title(tistr)
nexttile(TL) %B
figprefs.xlims = xLshort;
figprefs.doarrow = true;
figprefs.doleg = true; 
plotMdlDD(rpInfo,HLix,priColors,figprefs);

%%
h = .1; w = .1;
nR = TL.GridSize(1);
anns = cell(nR,1);
for r = 1:nR
    if r == 1; dy = 0.035; end 
    ann = annotation(gcf,'textbox',[0 1-((1/nR)*(r-1))-h-dy h h],...
        'Units','Normalized','String',['(',char('a'+(r-1)),')'],'EdgeColor','none','FontSize',10);
    ann.Position(1) = -0.015; % won't allow this in annotation arg
end

%%
testwrite = true;
if testwrite
    outdir = pwd;
	outfile = 'figS1_RP_no_oscil.tif';
	outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );
	fprintf('saving %s\n',strrep(outpath,'.tif','.eps'));
	print(gcf, strrep(outpath,'.tif','.eps'), '-depsc', '-r300' );  	
end %if