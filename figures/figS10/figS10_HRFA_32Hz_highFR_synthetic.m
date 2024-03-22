%figS10_HRFA_32Hz_highFR_synthetic.m
%%
clear all;  
close all;
clc;
warning('on','all')

scriptdirs = {pwd; '..\..\plot'; '..\..\analyze'; '..\..\helper_functions'};
addpath(scriptdirs{:})
 
datadir = pwd;
ssfile = 'subsample1000_highFR_stdRP.mat';   
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

fo = FOarr(end); % focus on 32 Hz

oti = sprintf('synthetic data: method hit rates and false alarms for f_{osc} = %u Hz (high FRs, %u ms RP, k = %0.1f)',fo,9,0.7);
otifs = 10;
nR = 4; nC = nNW;

%%
goArr = gobjects(nR*nC,1);
legArr = gobjects(nPBO,1);
k = 0;
cstrv = @(x) repmat(cellstr(x),nPBO,1);
modstr = 'modulation (%)';
XL.HR = [0 100]; 
YL.HR = [0 1];
XL.FA = [0 100];
YL.FA = [0 1]; 

gmin = +Inf;
gmax = -Inf;

%%
maxwidin = 7.5; maxhgtin = 8.75; hfac = 0.9; % 4 rows
figure('PaperPosition',[0 0 maxwidin maxhgtin*hfac])

%%
TL = tlcompact(nR,nC,oti,otifs);
codeArr = 'RC';
suffArr = {'C','R','DRC'};  
longSuffArr = {'comp','res','diff'};
ylblMap = containers.Map();
ylblMap('R') = 'res';
ylblMap('C') = 'shuf';
createNegated = true;
if createNegated
    dSuff = suffArr{end};
    for Var = {'HR','FA'}
        T.([Var{:},dSuff]) = -1*T.([Var{:},dSuff(1),dSuff(3),dSuff(2)]);
    end %for
end

%%
tvars = {'C','NW','FO','PBO','PB','M','HR'};
nV = numel(tvars);
vartypes = [{'char'};repmat({'double'},nV-1,1)];
nRtbl = numel(codeArr)*nNW*nPBO*nM; 
mnTbl = table('Size',[nRtbl nV],'VariableNames',tvars,'VariableTypes',vartypes);
rep6 = @(x) repmat(x,nM,1);

%%
dArr = {'HR','FA'};
rw = 1;

for c = 1:numel(codeArr)
    code = codeArr(c);
    for D = 1:numel(dArr)
        d = dArr{D};
        tfld = [d,code];
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
                uM = unique(m);
                mns = accumarray(ind,dvar,[],@mean);
                if isequal(d,'FA') 
                   tblrows = [rep6(nw),rep6(fo),rep6(pbo),rep6(pb),uM,mns];
                   mnTbl.C(rw:rw+nM-1) = rep6({code});
                   mnTbl(rw:rw+nM-1,2:end) = num2cell(tblrows);
                   rw = rw+nM;                
				end
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
                ax = gca;
                ax.XTick = [0:20:100];
                ax.XTickLabelRotation = 0;    
                ylim(yL) 
                ty = 0.93*yD+yL(1); dy = -diff(yL)/10; lx = 0.82*100;  
                lxd = lx-0.03*100;
                y = ty+p*dy;
                fs = 6;
                if n == 3 && ax.Layout.Tile <= 3 
                    if p == 1 
                        str = 'base FR (offset)';
                        ttxt = text(lxd,y-dy,str);
                        ttxt.FontWeight = 'bold';
                        ttxt.FontSize = fs;
                    end
                    str = sprintf('%u Hz (%u)',pb,pbo);
                    txt = text(lx,y,str);
                    txt.FontSize = fs;
                    txtd = text(lxd,y+0.01,'\cdot');
                    txtd.FontSize = fs+8;
                    txtd.Color = PL.Color;
                    txtd.FontWeight = 'bold';
                    txtd.VerticalAlignment = 'middle';
                end
                ax = gca;
                goArr(k) = ax;      
                legArr(p) = PL2;
                ax.TickDir = 'out';
                ax.Box = 'off';            
                if ax.Layout.Tile <= 3 
                    tstr = sprintf('duration = %u x 1024 ms',nw);
                    title(tstr)
                elseif ax.Layout.Tile == nR*nC - floor(nC/2) 
                    xlabel(modstr) 
                end %if
                if n == 1 
                    ylbl = sprintf('mean(%s(%s)) +/- SE',d,ylblMap(code));
                    ylabel(ylbl)
                end %if n == 1
            end %p    
        end %nw
    end %D
end %c
%%

%%		
testwrite = true;
if testwrite
    outdir = pwd;
    outfile = 'figS10_HRFA_32Hz_highFR_synthetic.tif';
    mkdir(outdir)
    outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );
	fprintf('saving %s\n',strrep(outpath,'.tif','.eps'));
	print(gcf, strrep(outpath,'.tif','.eps'), '-depsc', '-r300' );  
end
