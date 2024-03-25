%figS4_HR_diffs_pbase_vary_synthetic.m
%%
clear all;  
close all;
clc;
warning('on','all')
%%

scriptdirs = {pwd; '..\..\plot'; '..\..\analyze'; '..\..\helper_functions'};
addpath(scriptdirs{:})

datadir = pwd;
ssfile = 'subsample1000_pbase_vary.mat';   
sspath = fullfile(datadir,ssfile);

%%
fprintf('loading %s\n',sspath)
load(sspath,'T');

%%
NWarr = unique(T.NW);
nNW = numel(NWarr);
FOarr = unique(T.FO);
nFO = numel(FOarr);
PBOarr = unique(T.PBO);
nPBO = numel(PBOarr);
Marr = unique(T.M);
nM = numel(Marr);
nSS = unique(T.SS);
PBarr = unique(T.PB);
nPB = numel(PBarr);
%%

maxwidin = 7.5; maxhgtin = 8.75;
figure('PaperPosition',[0 0 maxwidin maxhgtin-(maxhgtin-3/8)/2])
%% 
codeStrMap = containers.Map();
codeStrMap('C') = 'shuffling';
codeStrMap('R') = 'residuals';

%%
T.HRDRC = T.HRDCR*-1;
T.FADRC = T.FADCR*-1;
%%
code = 'RC';
strA = codeStrMap(code(1));
strB = codeStrMap(code(2));
d = 'HR';
k = 0.7;
nr = 9;
dmap = containers.Map();
dmap('HR') = 'hit rate';
dmap('FA') = 'false alarm';
oti = sprintf('synthetic data: %s differences (residuals - shuffling), p_{base} varied directly (%u ms RP, k = %0.1f)',...
    dmap(d),nr,k);
otifs = 10;
nR = nPB; nC = nNW;
%%
TL = tlcompact(nR,nC,oti,otifs);
%%
goArr = gobjects(nR*nC,1);
legArr = gobjects(nPBO,1);
k = 0;
cstrv = @(x) repmat(cellstr(x),nPBO,1);
modstr = 'modulation (%)';
XL.HR = [0 100];
YL.HR = [-0.75 0.75];
XL.FA = [0 100];
YL.FA = [-0.75 0.75]; 

gmin = +Inf;
gmax = -Inf;
%%
tfld = [d,'D',code];
%%
for b = 1:nPB
    pb = PBarr(b);
    for n = 1:nNW
        nw = NWarr(n);
        nexttile
        k = k+1;
        for p = 1:nPBO
            pbo = PBOarr(p); 
            fo = pb-pbo;
            ix = T.NW == nw & T.PB == pb & T.PBO == pbo;
            m = T.M(ix)*100;
            ind = m/mode(diff(m))+1;
            dvar = T.(tfld)(ix);
            uM = unique(m);
            mns = accumarray(ind,dvar,[],@mean);
            if max(mns) > gmax; gmax = max(mns); end
            if min(mns) < gmin; gmin = min(mns); end
            se = accumarray(ind,dvar,[],@stderr);
            PL = plot(uM,mns);
            xL = XL.(d); yL = YL.(d);
            yD = diff(yL);
            hold on; 
            PL2 = plot(uM,mns,'.','Color',PL.Color,'LineStyle','none');
            e = errorbar(uM,mns,se,'Color',PL.Color);
            xlim(xL)
            xlim(xL)
            ax = gca;
            ax.XTick = [0:20:100];
            ax.XTickLabelRotation = 0;
            ylim(yL)     % for dy below was 11
            ty = 0.97*yD+yL(1); dy = -diff(yL)/10; lx = 0.82*100;
            lxd = lx-0.03*100;
            y = ty+p*dy;
            fs = 6;
            if n == nNW
                if p == 1
                    str = 'f_{osc} (offset)';
                    ttxt = text(lxd,y-dy,str);
                    ttxt.FontWeight = 'bold';
                    ttxt.FontSize = fs;
                end
                str = sprintf('%u Hz (%u)',fo,pbo);
                txt = text(lx,y,str);
                txt.FontSize = fs;
                %txtd = text(lxd,y,'\cdot');
                txtd = text(lxd,y+0.01,'\cdot');
                txtd.FontSize = fs+8;
                txtd.Color = PL.Color;
                txtd.FontWeight = 'bold';
                txtd.VerticalAlignment = 'middle';
            end   %n
            ax = gca;
            goArr(k) = ax;      
            legArr(p) = PL2;
            ax.TickDir = 'out';
            ax.Box = 'off';
            if b == 1
                tstr = sprintf('duration = %u x 1024 ms',nw);
                title(tstr)
            elseif b == nPB && n == ceil(nNW/2)
                xlabel(modstr)
            end %if
            if n == 1
                if b == ceil(nPB/2)
                    useFixedYlbl = true;
                    if ~useFixedYlbl
                        ylbl = sprintf('mean(HR(%s)-HR(%s)) +/- SE (note wider y limits)',strA,strB);
                    else
                        ylbl = sprintf('mean(%s(res)-%s(shuf)) +/- SE',d,d);
                    end
                    ylabel(ylbl)
                    if mod(nR,2) == 0
                        ax.YLabel.Position = [ax.YLabel.Position(1)+0.5 -0.85 ax.YLabel.Position(3)];
                    end
                end %f
                tx = text(0.05*100,0.75*yD+yL(1),sprintf('p_{base} = %u Hz',pb),'FontSize',7); 
            end %if
        end %p
    end %n
end %pb    

%%
testwrite = true;
if testwrite
    outdir = pwd;
    outfile = 'figS4_HR_diffs_pbase_vary_synthetic.tif';
    mkdir(outdir)
    outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );
	fprintf('saving %s\n',strrep(outpath,'.tif','.eps'));
	print(gcf, strrep(outpath,'.tif','.eps'), '-depsc', '-r300' );  
end
