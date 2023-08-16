%figS3_FA_diffs_primary_synthetic.m
%%
clear all;  
close all;
clc;
warning('on','all')

scriptdirs = {pwd; '..\..\plot'; '..\..\analyze'; '..\..\helper_functions'};
addpath(scriptdirs{:})
 
datadir = pwd;
ssfile = 'subsample1000_stdFR_stdRP.mat';   
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

%%
maxwidin = 7.5; maxhgtin = 8.75; 
figure('PaperPosition',[0 0 maxwidin maxhgtin])
%% 
codeStrMap = containers.Map();
codeStrMap('C') = 'shuffling';
codeStrMap('R') = 'residuals';

%%
code = 'CR';
negate = true; % so can see fldB-fldA
if negate
    strA = codeStrMap(code(2));
    strB = codeStrMap(code(1));
else
    strA = codeStrMap(code(1));
    strB = codeStrMap(code(2));
end

d = 'FA';
k = 0.7;
nr = 9;
dmap = containers.Map();
dmap('HR') = 'hit rate';
dmap('FA') = 'false alarm';

oti = sprintf('synthetic data: %s differences (residuals - shuffling), low to moderate FRs (%u ms RP, k = %0.1f)',...
    dmap(d),nr,k);
	
otifs = 10;
nR = nFO; nC = nNW;
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


for f = 1:nFO
	fo = FOarr(f);
	for n = 1:nNW
        nw = NWarr(n);
        nexttile
        k = k+1;
        for p = 1:nPBO  
            pbo = PBOarr(p);
            pb = pbo+fo;
            ix = T.NW == nw & T.FO == fo & T.PBO == pbo;
            m = T.M(ix)*100;
            ind = m/mode(diff(m))+1;
            dvar = T.(tfld)(ix);
            if negate; dvar = -1*(dvar); end
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
            ylim(yL)   
            ty = 0.97*yD+yL(1); dy = -diff(yL)/10; lx = 0.82*100;
            lxd = lx-0.03*100;
            y = ty+p*dy;
            fs = 6;
            if n == 3
                nrow = 2; spcx = 0.25;
                lx = (0.075 + spcx*idiv(p-1,nrow))*100;
                lxd = lx-0.03*100;
                ty = 0.3*yD+yL(1);
                dy = -diff(yL)/10;
                y = ty+(mod(p-1,nrow)+1)*dy;
                if p == 1
                    str = 'base FR (offset)';
                    ttxt = text(lxd,y-dy,str);
                    ttxt.FontWeight = 'bold';
                    ttxt.FontSize = fs;
                end %p
                str = sprintf('%u Hz (%u)',pb,pbo);
                txt = text(lx,y,str);
                txt.FontSize = fs;
                txtd = text(lxd,y+0.01,'\cdot');
                txtd.FontSize = fs+8;
                txtd.Color = PL.Color;
                txtd.FontWeight = 'bold';
                txtd.VerticalAlignment = 'middle'; 						
			end %if n
            ax = gca;
            goArr(k) = ax;      
            legArr(p) = PL2;
            ax.TickDir = 'out';
            ax.Box = 'off';			
            if f == 1
                tstr = sprintf('duration = %u x 1024 ms',nw);
                title(tstr)
            elseif f == nFO && n == ceil(nNW/2)
                xlabel(modstr)
            end %if
            if n == 1
                if f == ceil(nFO/2)
                    useFixedYlbl = true;
                    if ~useFixedYlbl
                        ylbl = sprintf('mean(FA(%s)-FA(%s)) +/- SE (note wider y limits)',strA,strB);
                    else
                        ylbl = sprintf('mean(%s(res)-%s(shuf)) +/- SE',d,d);
                    end
                    ylabel(ylbl)
                end %f
                tx = text(0.05*100,0.75*yD+yL(1),sprintf('f_{osc} = %u Hz',fo),'FontSize',7);
            end %if				
		end %p
	end %n
end %f


%%		
testwrite = true;
if testwrite
    outdir = pwd;
    outfile = 'figS3_FA_diffs_primary_synthetic.tif';
    mkdir(outdir)
    outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );  
end
