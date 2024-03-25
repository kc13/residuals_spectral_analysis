%figS13_eval_synthetic_RPk0_4nr4.m
%
%%
clear all;
close all;
clc;

%% addpath if not added yet
% assumes script run within local directory
scriptdirs = {pwd;
'..\..\plot';
'..\..\helper_functions';
'..\..\analyze';
}; 
addpath(scriptdirs{:})
datadir = pwd;
nSS = 1000;
ssfile = 'subsample1000_stdFR_RPk0.4nr4.mat';   
sspath = fullfile(datadir,ssfile);
gridfile = 'HRFAtable_stdFR_RPk0.4nr4.mat';
gridpath = fullfile(datadir,gridfile);
rpfile = 'RPtable_stdFR_RPk0.4nr4.mat';
rppath = fullfile(datadir,rpfile);

%%
fprintf('loading %s\n',sspath)
S = load(sspath);
fprintf('loading %s\n',gridpath)
G = load(gridpath);
fprintf('loading %s\n',rppath)
R = load(rppath);

%% fixed to highlight fig 2-3 examples
nw = 60;
f_osc = 9;
tbl = G.T;
ix = tbl.NW == nw & tbl.FO == f_osc;

suffArr = {'C','R','DRC'}; 
longSuffArr = {'comp','res','diff'};
createNegated = true; % res - comp
if createNegated
    dSuff = suffArr{end};
    for Var = {'HR','FA'}
        tbl.([Var{:},dSuff]) = -1*tbl.([Var{:},dSuff(1),dSuff(3),dSuff(2)]);
    end %for
end
subtbl = tbl(ix,:);
nG = numel(suffArr); % num grids
nM = S.nM; nPB = S.nPB;

%%
maxwidin = 7.5; maxhgtin = 8.75; 
figure('PaperPosition',[0 0 maxwidin maxhgtin])
oti = 'synthetic spike trains: method evaluation (low to moderate FRs, 4 ms RP, k = 0.4)';
otisz = 10.5;
nP = 3; nC = 1; 
TL = tlcompact(nP,nC,oti,otisz);
TL.TileIndexing = 'columnmajor';
%% 
colToMat = @(x) reshape(x,nPB,nM)';  %PB,nM
tOut = @(x) set(x,'TickDir','out');
pbCol = subtbl.PB;
pbMat = colToMat(pbCol);
mCol = subtbl.M*100;  %to show as %
mMat = colToMat(mCol);
pboCol = subtbl.PBO;
pboMat = colToMat(pboCol);

%%
nexttile(TL)
TL2 = tiledlayout(TL,2,3,'Padding','compact','TileSpacing','compact'); 
TL2.Layout.Tile = 1;
TL2.TileIndexing = 'columnmajor';
tstr = 'results sample #1: duration = 60 x 1024 ms, 9 Hz oscillation, 100 simulation iterations per cell';
TL2.Title.String = tstr;
TL2.Title.FontSize = 9;
axis off; box off;
modstr = 'modulation (%)';
pbstr = 'base FR (Hz)';
pbostr = 'base FR-oscil. freq.';


for g = 1:nG
    figopt = struct;
    figopt.doXlbl = false;
    if g==1
        figopt.doYlbl = true;
        figopt.ylbl = pbstr;
    elseif g==3
        figopt.doYRlbl = true;
        figopt.yrlbl = pbostr;
        figopt.yrmat = pboMat;
    end
	figopt.yrmat = pboMat;
	suff = suffArr{g}; 
	longsuff = longSuffArr{g};
    hitCol = subtbl.(['HR',suff]);
    hitMat = colToMat(hitCol);
    faCol = subtbl.(['FA',suff]);
    faMat = colToMat(faCol);  
	nexttile(TL2)
	ax = gridPlot(mMat,pbMat,hitMat,'hr',longsuff,figopt);
    nexttile(TL2)
    figopt = struct;
    figopt.doXlbl = true;
    figopt.xlbl = modstr;
    figopt.yrmat = pboMat;
    if g==1
        figopt.doYlbl = true;
        figopt.ylbl = pbstr;
    elseif g==3
        figopt.doYRlbl = true;
        figopt.yrlbl = pbostr;
        figopt.yrmat = pboMat;
    end
    ax = gridPlot(mMat,pbMat,faMat,'fa',longsuff,figopt);
end %g
pbo = 2;
m = 0.6;
ix = tbl.PBO == pbo & tbl.M == m;
subtblB = tbl(ix,:);	
nexttile(TL)
TLB = tiledlayout(TL,2,3,'Padding','compact','TileSpacing','compact'); 
TLB.Layout.Tile = 2;
TLB.TileIndexing = 'columnmajor';
tstr = 'results sample #2: base FR-oscil. freq. = 2 Hz, moduation = 60%, 100 simulation iterations per cell';
TLB.Title.String = tstr;
TLB.Title.FontSize = 9;
axis off; box off;
fostr = 'oscil. freq. (Hz)';
nwstr = 'duration (x 1024 ms)';
nNW = S.nNW; nFO = S.nFO;
colToMatB = @(x) reshape(x,nFO,nNW);
nwCol = subtblB.NW;
nwMat = colToMatB(nwCol);
foCol = subtblB.FO;
foMat = colToMatB(foCol);
pbColB = subtblB.PB;
pbMatB = colToMatB(pbColB);
for g = 1:nG
    figopt = struct;
    figopt.doXlbl = false;
    if g==1
        figopt.doYlbl = true;
        figopt.ylbl = fostr;
    elseif g==3
        figopt.doYRlbl = true;
        figopt.yrlbl = pbstr;
        figopt.yrmat = pbMatB;    
    end 
	suff = suffArr{g};
    longsuff = longSuffArr{g};
    hitCol = subtblB.(['HR',suff]);
    hitMat = colToMatB(hitCol); 
    faCol = subtblB.(['FA',suff]);
    faMat = colToMatB(faCol);
    nexttile(TLB)
    figopt.yrmat = pbMatB;   
	ax = gridPlot(nwMat,foMat,hitMat,'hr',longsuff,figopt);
    nexttile(TLB)
    figopt = struct;
    figopt.doXlbl = true;
    figopt.xlbl = nwstr;
    if g==1
        figopt.doYlbl = true;
        figopt.ylbl = fostr;
    elseif g==3
        figopt.doYRlbl = true;
        figopt.yrlbl = pbstr;
        figopt.yrmat = pbMatB;       
    end     
    figopt.yrmat = pbMatB;   
    ax = gridPlot(nwMat,foMat,faMat,'fa',longsuff,figopt);
end %g
% C: RP
nexttile(TL);
axis off; box off;
TLC = tiledlayout(TL,1,2,'Padding','compact','TileSpacing','compact'); 
TLC.Layout.Tile = 3;
nexttile(TLC)
allRP = R.T.RP;
H = histogram(allRP,[1:max(allRP)+1],'orientation','horizontal');
H.FaceColor = [0.0039 0.1255 0.2353];
xlbl = 'count (of 540 cells x 100 iterations)';
xlabel(xlbl)
ylbl = 'estimated duration (ms)';
ylabel(ylbl)
tstr = "recovery period estimates, target = " + num2str(S.simOpt.nr) + " ms";
med = median(allRP);
IQR = iqr(allRP);
txtstr = {sprintf('median = %u ms',med); sprintf('iqr = %u ms',IQR)};
ax = gca;
tX = ax.XLim(1) + 0.60*diff(ax.XLim);
tY = ax.YLim(1) + 0.75*diff(ax.YLim);
tx = text(tX,tY,txtstr,'FontSize',8);
axis tight
tOut(gca);
set(gca,'Box','off')
title(tstr)
% these lines obtain information
% about the RP distribution
absE = abs(allRP-R.simOpt.nr);
E = allRP-R.simOpt.nr;
dmn = @(x) x-mean(x);
pAcc = 1-[nnz(E)/numel(E)];
[ctsE,edgesE] = histcounts(E,min(E):1:max(E)+1);
edgesEL = edgesE(1:end-1);
pE = ctsE/numel(E);
[ctsAE,edgesAE] = histcounts(absE,min(absE):1:max(absE)+1);
edgesAEL = edgesAE(1:end-1);
pAE = ctsAE/numel(absE);
dispRPInfo = true;
maxAE = 15-9;
if dispRPInfo  
    for k = 0:maxAE   
        fprintf('%0.2f%% \t<= %u ms\n',round(sum(pAE(edgesAEL <= k))*100,2),k); 
    end  
    for k = 0:maxAE-1
        fprintf('%0.2f%%, ',round(sum(pAE(edgesAEL <= k))*100,2)); 
    end 
    k = k+1;
    fprintf('%0.2f%%\n',round(sum(pAE(edgesAEL <= k))*100,2)); 
end %if
% RP plot
nexttile(TLC)
axis off; box off;
tstr = 'partial ROC curves: bootstrapped pAUC comparison';
TLD = tiledlayout(TLC,1,3,'Padding','compact','TileSpacing','compact');
TLD.Layout.Tile = 2;
TLD.Title.String = tstr;
TLD.Title.FontSize = 9;

%pROC
nexttile(TLD,[1 2])
fld1List = S.fld1arr;
nF1 = numel(fld1List);
faMatCropFA = S.faMatCropFA;
hrMatCropFA = S.hrMatCropFA;
hold all;

gArr = gobjects(2,1);
tmpClrs = {[0.3020 0.0784 0.5490];[0.0039 0.1255 0.2353]};
nanArr = gobjects(2,1);  % for formatting purposes
for f1 = 1:nF1
    fld1 = fld1List{f1};
    faData = faMatCropFA.(fld1); 
    hrData = hrMatCropFA.(fld1);
    hold all;
    p = cellfun(@(x,y) plot(x',y','Marker','.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]),faData,hrData,'UniformOutput',false);
    for i = 1:numel(p)
        p{i}.Color = tmpClrs{f1};
    end
    gArr(f1) = p{1};
	% placeholder nan line to enable 
    % thicker formatting of legend lines
    nanline = plot(faData{1}',nan(size(hrData{1}))',...
        'LineWidth',p{1}.LineWidth*3,'Color',p{1}.Color,...
        'Marker','.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
    nanArr(f1) = nanline;
end %for
hold off;
ax = gca;
xlabel('false alarm rate')
ylabel('hit rate')
axis tight  
ax.Box = 'off';
ax.TickDir = 'out';
tstr = {'pROC: 1000 paired subsamples'};
title(tstr)
legend(nanArr,{'shuffling','residuals'},'Box','off','Location','Southeast')
tX = ax.XLim(1) + 0.15*diff(ax.XLim);
tY = ax.YLim(1) + 0.40*diff(ax.YLim);
txtstr = {['criterion = '],['significance test threshold']};
tx = text(tX,tY,txtstr,'FontSize',8);
nexttile(TLD)
S.DiffsF.RC = -1*S.DiffsF.CR;
S.mnDiffF.RC = mean(S.DiffsF.RC);
S.seDiffF.RC = stderr(S.DiffsF.RC);
[S.hF.RC,S.pF.RC,S.ciF.RC,S.statsF.RC]...
    = ttest(S.pAUCf.res,S.pAUCf.comp);
b = bar(S.mnDiffF.RC);
b.FaceColor = [0.5 0.5 0.5];
hold on;
e = errorbar(S.mnDiffF.RC,diff(S.ciF.RC)/2); 
e.LineWidth = 0.5;  % default 0.5
e.CapSize = 10;  
mn = S.mnDiffF.RC;
se = S.seDiffF.RC;
% some code adapted from Matlab's buiit-in ttest fn:
alpha = 0.01;
df = S.statsF.RC.df;
ser = S.statsF.RC.sd/sqrt(numel(S.DiffsF.RC));
crit = tinv((1 - alpha / 2), df) .* ser;
dim = 1;
xmean = mn;
ci = cat(dim, xmean - crit, xmean + crit);
e = errorbar(S.mnDiffF.RC,diff(ci)/2);
tstr = {'pAUC difference'};
title(tstr)
ax.Box = 'off';
ax.TickDir = 'out';
pd = 0.001;
ylim([0 ci(2)+pd])
ax = gca;
ax.Box = 'off';
ax.TickDir = 'out';
legend({'mean';'99% CI'},'Box','off','Location','southoutside')
xticklabels({''})
yticks([0,0.024])
ylabel('pAUC(res)-pAUC(shuf)')
%%

% panel labels
h = .1; w = .1;
nR = TL.GridSize(1);
for r = 1:nR
    if r == 1; dy = 0.035; elseif r == 2; dy = 0.005; else; dy = 0; end
    annotation(gcf,'textbox',[0 1-((1/nR)*(r-1))-h-dy h h],...
        'Units','Normalized','String',['(',char('a'+(r-1)),')'],'EdgeColor','none','FontSize',10)
end
annotation(gcf,'textbox',[0.5 1-((1/nR)*(r-1))-h h h],...
        'Units','Normalized','String',['(',char('a'+(r)),')'],'EdgeColor','none','FontSize',10)
		
%%		
testwrite = false;
if testwrite
    outdir = pwd;
    outfile = 'figS13_eval_synthetic_RPk0.4nr4.tif';
    mkdir(outdir)
    outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );
	fprintf('saving %s',strrep(outpath,'.tif','.eps'));
	print(gcf, strrep(outpath,'.tif','.eps'), '-depsc', '-r300' );  
end

%% ratios of hit rates and FAs at p < 0.05 corrected level
hrMat = S.hrMns;
faMat = S.faMns;
%%
aIX = S.aIX;
alphaArr = S.alphaArr;
hrVec05 = structfun(@(x) x(:,aIX),hrMat,'UniformOutput', false);
faVec05 = structfun(@(x) x(:,aIX),faMat,'UniformOutput', false);

%%
mnHR05 = structfun(@(x) mean(x), hrVec05, 'UniformOutput', false);
mnFA05 = structfun(@(x) mean(x), faVec05, 'UniformOutput', false);

%% confirm sig diffs
[hHV,pHV,ciHV,statsHV] = ttest(hrVec05.res,hrVec05.comp);
[hFV,pFV,ciFV,statsFV] = ttest(faVec05.res,faVec05.comp);

%% summary statement
fprintf('at alpha_c = %0.2f, residuals mean HR = %0.2f%%, shuffling mean HR = %0.2f%%\n',...
    alphaArr(aIX),mnHR05.res*100,mnHR05.comp*100);
fprintf('at alpha_c = %0.2f, residuals mean FA = %0.2f%%, shuffling mean FA = %0.2f%%\n',...
    alphaArr(aIX),mnFA05.res*100,mnFA05.comp*100);
