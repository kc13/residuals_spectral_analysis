% fig3_res_illustrate.m
% generates figure 2 using stored values used for manuscript
%%
clear all;
close all;
clc;

%% addpath if not added yet
% assumes script run within local directory
scriptdirs = {'../../plot';'../../helper_functions'};
addpath(scriptdirs{:})
datadir = '../fig2';

datafile = 'figs_2_3_shuf_res_illustrate_data.mat';
datapath = fullfile(datadir,datafile);
%%
colorfile = 'priColors.mat';
colordir = '../../plot';
colorpath = fullfile(colordir,colorfile);

%%
fprintf('loading %s\n',datapath)
load(datapath)
fprintf('loading %s\n',colorpath)
load(colorpath)

%%
maxwidin = 7.5; maxhgtin = 8.75;
figure('PaperPosition',[0 0 maxwidin maxhgtin])

%%
oti = sprintf(['synthetic spike trains: %u Hz oscillation' ...
    ' (residuals method, axis bounds vary to aid visualization)'],L.simOpt.f_osc);
otisz = 11;
TL = tlcompact(4,2,oti,otisz);
TL.TileIndexing = 'columnmajor';

%% arrays of high/low FR info
rpArr = [H.res_results.res.rpInfoS;L.res_results.res.rpInfoS];
soArr = [H.simOpt;L.simOpt];
mArr = [H.m,L.m];
resArr = {H.res_results;L.res_results};  
siArr = [H.sim_inputs; L.sim_inputs];

%% plotting prefs
xLshort = [1 20];
xLmore = [1 200];
xlbl = sprintf('ISI: zoom to [%u, %u] ms',xLmore);
HLix = [1 5 L.res_results.res.rpInfoS.optLag]; 

%%
nR = numel(rpArr);
figprefs = struct();
spkdArr = cell(nR,1);
fitpadArr = cell(nR,1);
respadArr = cell(nR,1);

%%
for r = 1:nR
    rpInfo = rpArr(r);
    simOpt =soArr(r);
    m = mArr(r);
    nexttile(TL) %A
    figprefs.xlims = xLmore;
    figprefs.xlbl = xlbl;
    if r == 1; figprefs.doleg = false; else; figprefs.doleg = true; end
	plotExpToISI(rpInfo,HLix,priColors,figprefs);
    tistr = sprintf('example #%u: base FR = %u Hz, modulation = %u%%',r,simOpt.pbase*1000,m*100);
	title(tistr)
	nexttile(TL) %B
    figprefs.xlims = xLshort;
    if r == 1; figprefs.doarrow = true; else; figprefs.doarrow = true; end
    if r == 1; figprefs.doleg = false; else; figprefs.doleg = true; end
	plotMdlDD(rpInfo,HLix,priColors,figprefs);
	nexttile(TL); %C
	ybds = [-.05 .1];
	tps = 'compact'; %'none'
    TL2 = tiledlayout(TL,1,3,'Padding','compact','TileSpacing',tps);
    TL2.Layout.Tile = 3+(r-1)*4; 
    axis off; box off;
	nexttile(TL2) % orig
    msecs = 1:209;
    si = siArr(r);
    res = resArr{r}.res;
    mdl = res.mdls1;
    spkd = si.src_spk_delta;
    spkdArr{r} = spkd;
	plotSpikeSample(spkd,msecs,'orig');
    ylim(ybds)  
    nexttile(TL2) % fit
	pd = rpInfo.RPend;
    fitpad = vertcat(nan(pd,1),mdl.Fitted.Response);
    fitpadArr{r} = fitpad;
    plotSpikeSample(fitpad,msecs,'fit');
    ylim(ybds)
    yticks('')
	nexttile(TL2) % res
    respad = vertcat(nan(pd,1),mdl.Residuals.Raw);
    respadArr{r} = respad;
    plotSpikeSample(respad,msecs,'res');
    ylim(ybds)
    yticks('')
    hold off;
	nexttile(TL,4+(r-1)*4) %D	
    tistr = sprintf('residuals-corrected (%u Hz/%u%%)',simOpt.pbase*1000,H.m*100);
    title(tistr)
    xticks(''); yticks(''); box off;
    TL3 = tiledlayout(TL,1,4,'Padding','compact','TileSpacing','none');
    TL3.Layout.Tile = 4+(r-1)*4;
    if r == 1; figopt.doleg = false; else; figopt.doleg = true; figopt.legloc = 'Southeast'; end
    figopt.plotShuf = false;
    figopt.arrowfreq = H.simOpt.f_osc; 
    figopt.inputax = nexttile(TL3);
    res = resArr{r};
    gArr = gobjects(2,1);
    splitF = 25;
    figopt.xlim = [0 splitF];
    figopt.XT = [0:5:splitF];
    figopt.doleg = false;
    axZ = specStatFig(res,'res','S1',figopt);
    gArr(1) = axZ;
    title('')
    xlabel('zoom 0-25 Hz')
    if r == 2
        figopt.doleg = true;
    end
    figopt.xlim = [splitF 500];
    figopt.inputax = nexttile(TL3,[1 3]);
    figopt.XT = [0:100:500];
    figopt.arrowfreq = [];
    ax = specStatFig(res,'res','S1',figopt);
    gArr(2) = ax;
    xlabel('frequency (25-500 Hz)')
    title('')
    yticks(''); ylabel('');
    matchYLims(gArr,1,[],[]);	
end % nR

%%
h = .1; w = .1;
nR = TL.GridSize(1);
for r = 1:nR
    annotation(gcf,'textbox',[0 1-((1/nR)*(r-1))-h h h],...
        'Units','Normalized','String',['(',char('a'+(r-1)),')'],'EdgeColor','none','FontSize',10)
end

%%
testwrite = true;
if testwrite
    outdir = pwd;
	outfile = 'fig3_res_illustrate.tif';
	outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );
	fprintf('saving %s\n',strrep(outpath,'.tif','.eps'));
	print(gcf, strrep(outpath,'.tif','.eps'), '-depsc', '-r300' );  	
end %if
