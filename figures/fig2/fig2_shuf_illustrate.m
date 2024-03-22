% fig2_shuf_illustrate.m
% generates figure 2 using stored values used for manuscript
%%
clear all;
close all;
clc;

%% addpath if not added yet
% assumes script run within local directory
scriptdirs = {'../../plot';'../../helper_functions'};
addpath(scriptdirs{:})
datadir = pwd;
datafile = 'figs_2_3_shuf_res_illustrate_data.mat';
datapath = fullfile(datadir,datafile);

%%
fprintf('loading %s\n',datapath)
load(datapath)

%%
maxwidin = 7.5; maxhgtin = 8.75;
figure('PaperPosition',[0 0 maxwidin maxhgtin])

%%
goArr = gobjects(2,1);
oti = sprintf(['synthetic spike trains: %u Hz oscillation' ...
    ' (shuffling method, axis bounds vary to aid visualization)'],L.simOpt.f_osc);
otisz = 11;
TL = tlcompact(4,2,oti,11); 
nexttile(TL,1) % panel A, High FR
tistr = sprintf('example #1: base FR = %u Hz, modulation = %u%%',H.simOpt.pbase*1000,H.m*100);
T = title(tistr);
tshift = [0.055 0 0];
T.Position = T.Position + tshift;	
xticks(''); yticks(''); box off;
TLH = tiledlayout(TL,1,4,'Padding','compact','TileSpacing','none');
TLH.Layout.Tile = 1;
figopt.doleg = false;
figopt.plotShuf = true;
figopt.arrowfreq = H.simOpt.f_osc;
figopt.inputax = nexttile(TLH);
splitF = 25;
figopt.xlim = [0 splitF];
figopt.XT = [0:5:25];
gArr = gobjects(2,1);
axHZ = specStatFig(H.shuf_results,[],'S1',figopt);
gArr(1) = axHZ;
title('')
xlabel('zoom 0-25 Hz')
goArr(1) = axHZ;
figopt.inputax = nexttile(TLH,[1 3]);
figopt.xlim = [splitF 500];
figopt.XT = [100:100:500];
figopt.arrowfreq = [];
axH = specStatFig(H.shuf_results,[],'S1',figopt);
xlabel('frequency (25-500 Hz)')
gArr(2) = axH;
matchYLims(gArr,1,[],[]);
title('')
yticks(''); ylabel('');
nexttile(TL,2) % panel A, Low FR
tistr = sprintf('example #2: base FR = %u Hz, modulation = %u%%',L.simOpt.pbase*1000,L.m*100);
T = title(tistr);
T.Position = T.Position + tshift;
xticks(''); yticks(''); box off;
TLL = tiledlayout(TL,1,4,'Padding','compact','TileSpacing','none');
TLL.Layout.Tile = 2;
figopt.doleg = false;
figopt.plotShuf = true;
figopt.arrowfreq = L.simOpt.f_osc;  % new 8/18/22
figopt.inputax = nexttile(TLL);
splitF = 25;
figopt.xlim = [0 splitF];
figopt.XT = [0:5:splitF];
axLZ = specStatFig(L.shuf_results,[],'S1',figopt);
gArr(1) = axLZ;
title('')
xlabel('zoom 0-25 Hz')
figopt.doleg = true;
figopt.legloc = 'Southeast';
figopt.xlim = [splitF 500];
figopt.inputax = nexttile(TLL,[1 3]);
figopt.XT = [100:100:500];
figopt.arrowfreq = [];
axL = specStatFig(L.shuf_results,[],'S1',figopt);
xlabel('frequency (25-500) Hz')
gArr(2) = axL;
matchYLims(gArr,1,[],[]);
title('')
yticks(''); ylabel('');
figopt = rmfield(figopt,{'xlim'});
% start here prep for panel B  
Lspkd = L.sim_inputs.src_spk_delta;
eIX = (L.simOpt.len-999:L.simOpt.len);
LsdE = Lspkd(eIX);
% 'Lsvec','Ldshuff','LrpIx','Lisi' loaded from mat
% generated during an earlier randomization of the demo ISIs
TL.Padding = 'tight';
nexttile(TL,3,[1 2])  
LstE = find(LsdE);
t = colormap(turbo);
nI = numel(LstE)-1;
tshort = t(flipdim(round(256:-size(t,1)/nI:1),2),:); % 
cmat = tshort;
nT = numel(LstE);
y = 0.05;
lw = 2;
x1 = 1;
x2 = LstE(1);
vclr = [0 0 0];
vmin = y;
vh = 0.05;
ty = y+0.01;
lbly = vmin+vh;
lblx = 2;
ftsz = 8; % defaults to 10 on some platforms
line([x1 x2],[y y],'Color',[0.5 0.5 0.5],'LineWidth',lw)
hold all;
for I = 1:nT-1
    x1 = LstE(I);
    x2 = LstE(I+1);
    line([x1 x2],[y y],'Color',cmat(I,:),'LineWidth',lw)
    text(mean([x1,x2]),y,num2str(Lisi(I)),'FontSize',ftsz,...
        'HorizontalAlignment','center','VerticalAlignment','top')  % was ty
    line([x1 x1],[vmin vmin+vh],'Color',vclr,'LineWidth',lw)
    %line([x1 x1],[vmin 0.2],'Color',cmat(I,:),'LineWidth',lw)
end %I
x1 = LstE(nT);
x2 = numel(LsdE);
line([x1 x2],[y y],'Color',[0.5 0.5 0.5],'LineWidth',lw)
line([x1 x1],[vmin vmin+vh],'Color',vclr,'LineWidth',lw)
grid on;
text(lblx,lbly,'original')
ylim([-0.1 0.15])
%% next layer, stacked below
st = find(Lsvec);
y2 = -0.05;
vmin2 = y2;
ty2 = y2+0.01;
x1 = 1;
x2 = st(1);
lbly = vmin2+vh;
line([x1 x2],[y2 y2],'Color',[0.5 0.5 0.5],'LineWidth',lw)
for I = 1:nT-1
    x1 = st(I);
    x2 = st(I+1);
    line([x1 x2],[y2 y2],'Color',cmat(LrpIx(I),:),'LineWidth',lw)
    text(mean([x1,x2]),y2,num2str(Ldshuff(I)),'FontSize',ftsz, ...
        'HorizontalAlignment','center','VerticalAlignment','top') % was ty2
    line([x1 x1],[vmin2 vmin2+vh],'Color',vclr,'LineWidth',lw)
end %I
x1 = st(nT);
x2 = numel(Lsvec);
line([x1 x2],[y2 y2],'Color',[0.5 0.5 0.5],'LineWidth',lw)
line([x1 x1],[vmin2 vmin2+vh],'Color',vclr,'LineWidth',lw)
grid on;
text(lblx,lbly,'shuffled')
ylim([-0.1 0.15])
yticks('')
xlabel('msec')
set(gca,'XColor','none','YColor','none')
grid off
title('illustration of ISI shuffling')
%% panel C shuffling (comp)ensation
nexttile(TL,5) % panel C, High FR
Htistr = sprintf('shuffle-corrected (%u Hz/%u%%)',H.simOpt.pbase*1000,H.m*100);
title(Htistr)
xticks(''); yticks(''); box off;
TLH = tiledlayout(TL,1,4,'Padding','compact','TileSpacing','none');
TLH.Layout.Tile = 5;
figopt.doleg = false;
figopt.plotShuf = false;
figopt.arrowfreq = H.simOpt.f_osc;  % new 8/18/22
figopt.inputax = nexttile(TLH);
splitF = 25;
figopt.xlim = [0 splitF];
figopt.XT = [0:5:25];
gArr = gobjects(2,1);
axHZ = specStatFig(H.shuf_results,'comp','S1Comp',figopt); % left zoom
gArr(1) = axHZ;
title('')
xlabel('zoom 0-25 Hz')
goArr(1) = axHZ;
figopt.inputax = nexttile(TLH,[1 3]);
figopt.xlim = [splitF 500];
figopt.XT = [0:100:500];
xlabel('25-500 Hz')
figopt.arrowfreq = [];
axH = specStatFig(H.shuf_results,'comp','S1Comp',figopt);
gArr(2) = axH;
matchYLims(gArr,1,[],[])
axH = specStatFig(H.shuf_results,[],'S1',figopt);
title('')
xlabel('frequency (25-500 Hz)')
yticks(''); ylabel('');
nexttile(TL,6) % panel C, L
Ltistr = sprintf('shuffle-corrected (%u Hz/%u%%)',L.simOpt.pbase*1000,L.m*100);
title(Ltistr)
xticks(''); yticks(''); box off;
TLL = tiledlayout(TL,1,4,'Padding','compact','TileSpacing','none');
TLL.Layout.Tile = 6;
figopt.doleg = false;
figopt.plotShuf = false;
figopt.arrowfreq = L.simOpt.f_osc;  
figopt.inputax = nexttile(TLL);
figopt.xlim = [0 splitF];
figopt.XT = [0:5:25];
axLZ = specStatFig(L.shuf_results,'comp','S1Comp',figopt);
gArr(1) = axLZ;
title('')
xlabel('zoom 0-25 Hz')
goArr(2) = axLZ;
figopt.doleg = true;
figopt.legloc = 'Southeast';
figopt.xlim = [splitF 500];
figopt.XT = [100:100:500];
figopt.inputax = nexttile(TLL,[1 3]);
figopt.arrowfreq = [];
axL = specStatFig(L.shuf_results,'comp','S1Comp',figopt);
xlabel('frequency (25-500 Hz)')
% little lower
Leg = axL.Parent.Children(1);
Leg.Position = Leg.Position + [0.2 -0.25 0 0];
gArr(2) = axL;
matchYLims(gArr,1,[],[])
title('')
yticks(''); ylabel('');
%% panel D : ISI PDFs
% using res result stucture because it stores this ISI info
goArr = gobjects(2,1);
rpiL = L.res_results.res.rpInfoS;
rpiH = H.res_results.res.rpInfoS;
pdfL = rpiL.isiPDF;
pdfH = rpiH.isiPDF;
lagsL = rpiL.isiLags;
lagsH = rpiH.isiLags;
sp = 0.15; 
nexttile(TL,7)
bar(lagsH,pdfH,'FaceColor','k');
hold on;
f = fit(lagsH,pdfH,'smoothingspline','SmoothingParam',sp);
plot(f,lagsH,pdfH)
ax = gca;
delete(ax.Children(2))
lnclr = [249 99 2]/255; 
lnw = 0.6;
ax.Children(1).Color = lnclr;
ax.Children(1).LineWidth = lnw;
Htistr = sprintf('ISI distribution (%u Hz/%u%%)',H.simOpt.pbase*1000,H.m*100);
title(Htistr)
legend off
goArr(1) = gca;
nexttile(TL,8)
bar(lagsL,pdfL,'FaceColor','k');
hold on;
f = fit(lagsL,pdfL,'smoothingspline','SmoothingParam',sp);
p = plot(f,lagsL,pdfL);
ax = gca;
delete(ax.Children(2))
ax.Children(1).Color = lnclr;
ax.Children(1).LineWidth = lnw;
Ltistr = sprintf('ISI distribution (%u Hz/%u%%)',L.simOpt.pbase*1000,L.m*100);
title(Ltistr)
legend(p(2),'smoothed','Location','Northeast','Box','off')
goArr(2) = gca;
ax = gca;
for g = 1:numel(goArr)
    ax = goArr(g);
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.XLabel.String = 'lags (ms)';
    ax.YLabel.String = 'p(ISI = lag)';
    ax.YLim(1) = 0;
end
%% add panel labels
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
    outfile = 'fig2_shuf_illustrate.tif';
    outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );
	fprintf('saving %s\n',strrep(outpath,'.tif','.eps'));
	print(gcf, strrep(outpath,'.tif','.eps'), '-depsc', '-r300' );  
end
