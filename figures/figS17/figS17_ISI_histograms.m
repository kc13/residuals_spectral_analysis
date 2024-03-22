%fig_S17_ISI_histograms.m

%%
clear all;  
close all;
clc;
warning('on','all')

%% addpath if not added yet
% assumes script run within local directory
scriptdirs = {pwd;
'..\..\plot';
'..\..\helper_functions';
'..\..\analyze'
}; 
addpath(scriptdirs{:})
%%
datadir = pwd;
datafile = 'fig5_figS16_unit_ex.mat'; 
datapath = fullfile(datadir,datafile);

%%
fprintf('loading %s\n',datapath)
load(datapath);

%%
deltadir = pwd;
deltafile = 'delta_mats_GPi_VLa.mat'; 
deltapath = fullfile(deltadir,deltafile);

%%
fprintf('loading %s\n',datapath)
load(datapath);
%% new 6/9
fprintf('loading %s\n',deltapath)
D = load(deltapath);
%%
colordir = '../../plot';
colorfile = 'priColors.mat';
colorpath = fullfile(colordir,colorfile);

%%
load(colorpath);

%%
maxwidin = 7.5; maxhgtin = 8.75; hfac = 0.85;
figure('PaperPosition',[0 0 maxwidin maxhgtin*hfac])

%%
oti = 'post-MPTP NHP spike trains: ISI histograms and RP estimates (axis bounds may vary)';
otisz = 10.5; 
nR = 3; nC = 2;
TL = tlcompact(nR,nC,oti,otisz);
TL.TileIndexing = 'columnmajor';

%%
rArr = [E(3).res_results; E(4).res_results]; 
rpArr = [rArr(1).res.rpInfoS; rArr(2).res.rpInfoS];
tiArr = {['example #1: ', aC{:}, ' unit'],...
    ['example #2: ', aD{:}, ' unit']};
uArr = [E(3).unit_data; E(4).unit_data];
sdArr = arrayfun(@(x) cell2mat(D.corr_results.(x.brain_area{:})(x.su_row).src_delta),uArr,'UniformOutput',false);

%% prefs
xLshort = [1 20];
xLmore = [1 200];
xL2D = [1 50];
xlbl = sprintf('p(ISI = lag): zoom to [%u, %u] ms',xLmore);
xlbl2D = sprintf('ISI(n): zoom to [%u, %u] ms',xL2D);
ylbl2D = sprintf('ISI(n+1)',xL2D);
clims = [0 0.0270]; % this was ix 27 auto clim
%%
priColorsNew = priColors;
clrsHL = priColors('HL');
priColorsNew('HL') = clrsHL(end,:);
nR = numel(rpArr);
figprefs = struct();
ax2D = gobjects(nR,1);

for r = 1:nR
    rpInfo = rpArr(r);
    nexttile(TL) 
    figprefs = struct;
    figprefs.xlims = xLmore;
    figprefs.xlbl = xlbl;
    figprefs.ylims = [0 max(rpInfo.isiPDF)*1.05];
    if r == 1; figprefs.doarrow = false; else; figprefs.doarrow = false; end
    if r == 1; figprefs.doleg = true; else; figprefs.doleg = true; end
    figprefs.doTri = true;
    plotOptToISI(rpInfo,priColors,figprefs);
    title(tiArr{r})
    ax = gca;
    nexttile
    HLix = rpInfo.optLag;
    figprefs.xlims = xLshort;
    plotMdlDD(rpInfo,HLix,priColors,figprefs);
    if r == 2
        ax = gca;
        ax.Parent.Children(1).Location = 'northeast';
    end
    % 2D ISI
    udata = uArr(r);
    ix = uArr.su_row;
    sd = sdArr{r};
    [pdf2D,lags2D] = isiPDF2D(sd);
    nexttile(TL)
    figprefs.xlims = xL2D;
    figprefs.ylims = xL2D;
    figprefs.xlbl = xlbl2D;
    figprefs.ylbl = ylbl2D;
    figprefs.clims = clims;
    ax2D(r) = plotPDF2D(lags2D,lags2D,pdf2D,figprefs); 
end %r
%%
%% panel labels
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
    outfile = 'figS17_ISI_histograms.tif';
    mkdir(outdir)
    outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );
	fprintf('saving %s\n',strrep(outpath,'.tif','.eps'));
	print(gcf, strrep(outpath,'.tif','.eps'), '-depsc', '-r300' );  
end

%%
