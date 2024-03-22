%figS11_highFR_FA_examples.m
%%
clear all;  
close all;
clc;
warning('on','all')

%%
scriptdirs = {pwd;
 '..\..\plot';
 '..\..\helper_functions';
 '..\..\analyze'};
addpath(scriptdirs{:})

%%
datadir = pwd;
LFfile = 'len122880f_osc7iter1.mat';
LFpath = fullfile(datadir,LFfile);
HFfile = 'len122880f_osc32iter1.mat';
HFpath = fullfile(datadir,HFfile);
%%
fprintf('loading %s\n',LFpath)
L = load(LFpath);
fprintf('loading %s\n',HFpath)
H = load(HFpath);

%%
pbArr = [92,117];
Lstruct = L.sim_data.(['pb',num2str(pbArr(1))]).m100;
Hstruct = H.sim_data.(['pb', num2str(pbArr(2))]).m100;
siArr = [Lstruct.sim_inputs_all{:}; Hstruct.sim_inputs_all{:}];

%%
stArr = {Lstruct,Hstruct};
tpfxArr = {'residuals','shuffling'};
foArr = [7,32];
f1arr = {'res','comp'};
f2arr = {'S1','S1Comp'};
fvec = [0:500/512:500]';
XL = [0 500];
XTdt = 50;
flbl = 'frequency (Hz)';
colArr = [1,1];
srch_bnds = [0.1 100];
cntl_bnds = [250 500]; 
Alpha = 0.05;

%%
maxwidin = 7.5; maxhgtin = 8.75;
figure('PaperPosition',[0 0 maxwidin maxhgtin])  

%%
oti = sprintf('synthetic spike trains: false alarm examples, high base FRs (m = %u%%, T = %u x 1024 ms)',100,120);
TL = tlcompact(16,2,oti,11);
TL.TileIndexing = 'columnmajor';

%%
for s = 1:numel(stArr)
    nexttile([4 1])
    st = stArr{s};
    f_osc = foArr(s);
    out5Ix = find(abs(fvec-f_osc)>5);
    f1 = f1arr{s};
    fld1 = f1;
    f2 = f2arr{s}; 
    f = st.f;
    figopt = struct;
    figopt.XL = XL;
    figopt.XTdt = XTdt;
    figopot.flbl = flbl;
    figopt.arrowfreq = f_osc; 
    specStruct = struct();
    specStruct.f = st.f;
    c = colArr(s);
    specStruct.(fld1).(f2) = st.(fld1).S1(:,c);
    specStruct.(fld1).stats.(f2) = get_spec_stats(specStruct.(fld1).(f2),specStruct.f,srch_bnds,cntl_bnds,Alpha);
    ax = specStatFig(specStruct,fld1,f2,figopt);
    tstr = sprintf('%s (p_{base} = %u Hz, f_{osc} = %u Hz)',tpfxArr{s},pbArr(s),f_osc);
	title(tstr,'FontSize',8.8)
    ax.YLim = [0 max(ax.Children(end).YData)*1.01];
    ax.Children(1).Position = ax.Children(1).Position - [0 0.05 0];
    si = siArr(s);	
    spkd = si.src_spk_delta;
    spr = si.src_sim_p_raw;
    sp = si.src_sim_p;
    rp = sp./spr;
    msecs = 1:209; % for now
    nexttile([3 1])
    ax1 = plotSpikeModelSample(spr,msecs,'oscil');
    nexttile([3 1])
    ax2 = plotSpikeModelSample(rp,msecs,'rp');
    nexttile([3 1])
    ax3 = plotSpikeModelSample(sp,msecs,'pspk');
    nexttile([3 1])
    ax4 = plotSpikeModelSample(spkd,msecs,'spikes');	
end %s

%% panel labels
h = .1; w = .1;
nR = 4; 
for r = 1:2
    annotation(gcf,'textbox',[0 1-((1/nR)*(r-1))-h h h],...
        'Units','Normalized','String',['(',char('a'+(r-1)),')'],'EdgeColor','none','FontSize',10)
end

%%		
testwrite = true;
if testwrite
    outdir = pwd;
    outfile = 'figS11_highFR_FA_examples.tif';
    mkdir(outdir)
    outpath = fullfile(outdir,outfile);
    set(gcf,'color','w');
    set(gcf,'InvertHardCopy','off')
    fprintf('saving %s\n',outpath)
    print(gcf, outpath, '-dtiff', '-r300' );
	fprintf('saving %s\n',strrep(outpath,'.tif','.eps'));
	print(gcf, strrep(outpath,'.tif','.eps'), '-depsc', '-r300' );  
end
