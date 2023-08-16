%figS16_additional_unit_ex.m
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
maxwidin = 7.5; maxhgtin = 8.75;    
figure('PaperPosition',[0 0 maxwidin maxhgtin])
%%
oti = ['post-MPTP NHP spike trains: candidate \beta oscillations ' ...
    '(axis bounds vary to aid visualization)'];
otisz = 10.5;
nR = 2; nC = 1; 
TL = tlcompact(nR,nC,oti,otisz);
TL.TileIndexing = 'rowmajor'; 
%%
xcArr = [E(3).corr_results; E(4).corr_results];
sArr = [E(3).shuf_results; E(4).shuf_results];
rArr = [E(3).res_results; E(4).res_results];
uArr = [E(3).unit_data; E(4).unit_data];
leglocArr = {'Southeast','Northeast'};
tiArr = {['example #1: ', aC{:}, ' unit'],...
    ['example #2: ', aD{:}, ' unit']};

%%
nE = numel(xcArr);
for e = 1:nE 
    NT = nexttile(TL);
    axis off; box off;
    quadopt = struct;
    quadopt.oti = tiArr{e};
    quadopt.otifs = 9;
    quadopt.XL = [0 100]; 
    quadopt.XTdt = 10;
    quadopt.flbl = 'frequency (Hz): zoom to [0, 100] Hz';
    quadopt.parent = TL; 
    quadopt.tilenum = e; 
    quadopt.plotShuf = true;
    quadopt.legloc = leglocArr(e);
    xcData = xcArr(e); sData = sArr(e); rData = rArr(e);
    fig = specStatQuad(xcData,sData,rData,quadopt); 
    if e == 2
        ax = gca;
        L = ax.Parent.Children(3);
        L.Visible = 'off';
    end
end %e
%%
F = gcf;
h = .1;
dy = [0, 0.01];
for r = [1:nR]
    r2 = nR-r+1;
    pos = F.Children.Children(r*2-1).Position;
    annotation(F,'textbox',[0 pos(2)+pos(4)+dy(r) h h],...
        'Units','Normalized','String',char('A'+(r2-1)),...
        'EdgeColor','none','FontSize',10)
end
%%
testwrite = true;
if testwrite
    outdir = pwd;
    outfile = 'figS16_additional_unit_ex.tif';
    mkdir(outdir)
    outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );  
end
%%


