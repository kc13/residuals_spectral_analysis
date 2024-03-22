%fig6_unit_population.m
% this will both recreate figure 6
% and run some descriptive and
% inferential stats 
% (see descriptions towards the end of the file)

%%
clear all;  
close all;
clc;
warning('on','all')

%% addpath if not added yet
% assumes script run within local directory
scriptdirs = {pwd;
'..\..\plot';
'..\..\analyze';
}; 
addpath(scriptdirs{:})

%%
datadir = pwd;
datafile = 'fig6_unit_population.mat'; 
datapath = fullfile(datadir,datafile);

%%
fprintf('loading %s\n',datapath)
load(datapath);
%%
maxwidin = 7.5; maxhgtin = 8.75;   
fracX = 1; fracY = 0.75; % for now, can't decide...
figure('PaperPosition',[0 0 maxwidin*fracX maxhgtin*fracY])
oti = 'post-MPTP NHP spike trains: GPi and VLa, thresholded PSDs (p_{corr} < 0.05)';
otisz = 10.5;  %this font size is from fig 3
TL = tlcompact(1,3,oti,otisz);
nF = numel(ff);
ln = -1*ones(1,nF);
flbl = 'frequency (Hz): zoom to [0, 100] Hz';
tstrArr = {'uncorrected','shuffling','residuals'};

%%
aCts = structfun(@(x) numel(x),FR.res);
getIM = @(x) labelPeaks(x,2) + x;
lnclr = ones(1,3)*0.5;
zclr = [1 1 1];
sclr = [0 0 1];  %sig ~pk
pclr = [1 0 1];  %peak
cmap = [lnclr;zclr;sclr;pclr];

%%
for f1 = 1:numel(f1strs)
	fld = f1strs{f1};
    [sGPi,ixGPi] = sort(FR.(fld).GPi); [sVLa,ixVLa] = sort(FR.(fld).VLa);
    m = vertcat(getIM(sigix.(fld).GPi(ixGPi,:)),ln,getIM(sigix.(fld).VLa(ixVLa,:)));
    frVec = vertcat(sGPi,nan,sVLa);
    nexttile
    colormap(cmap)	
	im = image(m+2,'CDataMapping','direct');
    title(tstrArr(f1))
    if f1 == 2
        xlabel(flbl)
    end
    ax = gca;
    ax.XLim = [0 100];
    pos = ax.Position;
    ax.TickDir = 'out';
    xticks([0:20:100])
    xtickangle(0)	
    if f1 == 1
        axL = gca;
        yLL = axL.YLim;
        axL.YTick = [aCts(1)/2,aCts(1)+aCts(2)/2];
        axL.YTickLabel = {'GPi unit (ascending FR)','VLa unit (ascending FR)'};
        axL.TickLength = [0 0];
        ytickangle(90);
    else
        axL = gca;
        axL.TickLength = [0 0];
        yticklabels('')
    end		
end %f1
ax3 = gca;
pos = ax3.Position;
fw = 'bold';
an = annotation('textbox',[pos(1)+0.05,pos(2)-0.09,pos(3)/3.8,pos(4)/26],...
    'BackgroundColor','w','String','non-sig.','FontSize',8.5,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Margin',0,'FontWeight',fw);
gp = 0.01; 
x2 = an.Position(1)+an.Position(3)+gp;
y2 = an.Position(2); w2 = an.Position(3); h2 = an.Position(4);
an2 = annotation('textbox',[x2 y2 w2 h2],...
    'BackgroundColor',sclr,'String','sig.','FontSize',8.5,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Margin',0,'FontWeight',fw);
x3 = an2.Position(1)+an2.Position(3)+gp;
an2 = annotation('textbox',[x3 y2 w2 h2],...
    'BackgroundColor',pclr,'String','sig. peak','FontSize',8.5,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Margin',0,'FontWeight',fw);
ax = gca;
ax.FontSize = 9.5;

testwrite = true;
if testwrite
	outdir = pwd;
	outfile = 'fig6_unit_population.tif';
	mkdir(outdir)
	outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );
	fprintf('saving %s\n',strrep(outpath,'.tif','.eps'));
	print(gcf, strrep(outpath,'.tif','.eps'), '-depsc', '-r300' );  
end

%% start of analysis sections

%% firing rate descriptives by region
strs = structfun(@(x) sprintf('%0.2f, %0.2f, %0.2f',round(min(x),2),round(median(x),2),round(max(x),2)),FR.orig,'UniformOutput',false);
[p,h,stats] = ranksum(FR.orig.GPi,FR.orig.VLa);
disp('FR info (Hz)')
disp(strs)

%% duration info, 1 tr = 1s
smry = @(x) [min(x) median(x) max(x)];
disp('duration info')
nA = numel(area_list);
for a = 1:nA
    Area = area_list{a};
    fprintf('%s %u,%0.1f,%u \n',Area,smry(nTr.res.(Area)))
end

%% for the proportions of units with alpha-beta power
% (for descriptives and input to McNemar's test in R)
% the generated table also summarizes information 
% about the estimated RPs

% first building a table:
tvars = {'Area','rAnyAB','sAnyAB','FR','RP'};
vartypes = {'cell','logical','logical','double','double'};
nV = numel(tvars);
nU = structfun(@(x) size(x,1),sigix.orig);
nR = sum(nU);
T = table('Size',[nR nV],'VariableNames',tvars,'VariableTypes',vartypes);

k = 1;
for a = 1:nA
    Area = area_list{a};
	T.Area(k:k+nU(a)-1) = repmat({Area},nU(a),1);
	rsig = sigix.res.(Area);
	ssig = sigix.comp.(Area);
	abIX = (ff >= 8 & ff <= 30)';
	rAB = rsig(:,abIX);
	sAB = ssig(:,abIX);
	rAnyAB.(Area) = any(rAB,2);
	sAnyAB.(Area) = any(sAB,2);
	T.rAnyAB(k:k+nU(a)-1) = rAnyAB.(Area);
	T.sAnyAB(k:k+nU(a)-1) = sAnyAB.(Area);
	T.RP(k:k+nU(a)-1) = RP.res.(Area);
	T.FR(k:k+nU(a)-1) = FR.orig.(Area);
	k = k+nU(a);
end %a

%% now finding units for which at least one 
% method flagged alpha-beta points
ixOR = T.rAnyAB | T.sAnyAB; 
rSel = T.rAnyAB(ixOR);
sSel = T.sAnyAB(ixOR);
dSel = rSel-sSel;
rOnly = dSel > 0;
frSel = T.FR(ixOR); % 1 if residuals only, 0 otherwise

%% logistic regression: does FR predict rOnly status?
X = frSel;
Y = rOnly;
[b,dev,stats] = glmfit(X,Y,'binomial'); % logit is default link
fprintf('for the logistic regression: B1 = %0.4f, dfe = %u, p(B1) = %0.5g\n',...
b(2),stats.dfe,stats.p(2))

%% table that will be used for the McNemar test in R
% the table reports counts for unit alpha-beta presence
% for the following combinations:
% (res+,shuf+), (res+,shuf-)
% (res-,shuf+), (res-,shuf-)
mtbl = [T.rAnyAB';~T.rAnyAB']*[T.sAnyAB,~T.sAnyAB];