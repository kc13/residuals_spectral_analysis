% run_subsample_and_ROC_pbase_vary.m
% note: random seeds not available for this case
%%
clear all;
close all;
clc;
%%
%% add fn paths if not added yet
% assumes script run within local directory
scriptdirs = {
'../../figures';
'../../analyze';
'../../helper_functions';
};

datadir = pwd;
%%
addpath(scriptdirs{:})
addpath(datadir)

%% set write prefs down at bottom
outdir = pwd;
suffix = 'pbase_vary_new'; 

%% no random seed option (original random 
% seeds not located)

%% some other relevant vars will load with data files
simOpt.win_len_msec = 1024;
simOpt.nr = 9;
simOpt.k = 0.7;
m_arr = [0:0.2:1]; 
nM = numel(m_arr);
nIter = 100; 
veclen = (2^nextpow2(simOpt.win_len_msec))/2+1;
fld1arr = {'comp','res'}; 
nF1 = numel(fld1arr);
codes = ['C','R']';
codepairs = nchoosek(codes,2);
fldpairs = nchoosek(fld1arr,2);
nCO = size(codepairs,1);
pbase_offset = -1*[2:2:8];  % if applying to pbase
nPB = numel(pbase_offset);
nPBO = numel(pbase_offset);

%% outer loops for len, f_osc
nwArr = [30;60;120];
nNW = numel(nwArr);
pbase_arr = [15 20]/1000;
nPB = numel(pbase_arr); nFO = nPBO;

%% for stats
cntl_bnds = [250 500];
srch_bnds = [0 100];
alphaMat = [0.01 0.05]'*(10.^[-6:1:2]);
alphaArr = reshape(alphaMat,[],1);
alphaArr(end) = [];
nA = numel(alphaArr);
fvec = [0:500/512:500]';
f_ix = fvec > srch_bnds(1) & fvec <= srch_bnds(2);

%% for subsamples
nSS = 1000; 
ssSZ = 20;
ssIXmat = cell2mat(arrayfun(@(x) randsample(nIter,ssSZ)',[1:nSS]','UniformOutput',false));

%% table preallocate
A = {'SS','NW','FO','PBO','PB','M'};
B = [join([repmat({'HR'},nF1,1) cellstr(codes)],'')' join([repmat({'FA'},nF1,1) cellstr(codes)],'')'];
C = [join([repmat({'HRD'},nCO,1) cellstr(codepairs)],'')' join([repmat({'FAD'},nCO,1) cellstr(codepairs)],'')'];
tvars = [A,B,C];
nV = numel(tvars);
nV = numel(tvars);
vartypes = repmat({'double'},1,nV);
nR = nNW*nFO*nPB*nM*nSS;
T = table('Size',[nR nV],'VariableNames',tvars,'VariableTypes',vartypes);
targA = 0.05;
k = 0;
aIX = find(alphaArr == targA);
%%
dispmod = 100;
tic;
for NW = 1:nNW
    nw = nwArr(NW);
    simOpt.len = nw*simOpt.win_len_msec;   
    for PB = 1:nPB
        FN = (NW-1)*nPB + PB;
        disp(FN)
        simOpt.pbase = pbase_arr(PB);
        pbase = simOpt.pbase;  
        foArr = simOpt.pbase*1000 + pbase_offset;
        fname = sprintf('len%upb%uiter%u.mat',simOpt.len,simOpt.pbase*1000,nIter);
        fpath = fullfile(datadir,fname);
        fprintf('loading %s\n',fpath);
        load(fpath)
        for s = 1:nSS
            if mod(s,dispmod) == 0
                fprintf('running ss %d of %d\n',s,nSS)
            end
            ssIX = ssIXmat(s,:);
            for f1 = 1:nF1 
                fld1 = fld1arr{f1};
                hitMat.(fld1) = nan(nFO,nM,nA);
                faMat.(fld1) = nan(nFO,nM,nA);
            end %f1            
            for FO = 1:nFO
                simOpt.f_osc = foArr(FO);
                f_osc = simOpt.f_osc;
                targIx = knnsearch(fvec,f_osc,'K',3)';
                out5Ix = find(abs(fvec-f_osc)>5);
                pbSec = simOpt.pbase*1000;
                fldstrFO = sprintf('fo%g',f_osc);
                pbo = pbSec-f_osc;
                    for M = 1:nM
                        m = m_arr(M);
                        fldstrM = sprintf('m%g',m*100); 
                        freq = sim_data.(fldstrFO).(fldstrM).f;
                        for f1 = 1:nF1
                            fld1 = fld1arr{f1};
                            S1 = sim_data.(fldstrFO).(fldstrM).(fld1).S1;
                            sig = nan(numel(freq),ssSZ);
                            p_osc = sim_data.(fldstrFO).(fldstrM).p_osc;
                            for a = 1:nA
                                Alpha = alphaArr(a);
                                for si = 1:ssSZ  % loop over 20
                                    i = ssIX(si);
                                    pow = S1(:,i);
                                    stats = get_spec_stats(pow,freq,srch_bnds,cntl_bnds,Alpha);
                                    sig(:,si) = stats.sigvec1TL; % fixed 6/26/22 was i
                                end %si
                                sigrows = sig(targIx,:);
                                hits = ~isnan(sigrows);
                                anyHits = any(hits,1);
                                hitRate = sum(anyHits)/ssSZ;  
                                if m == 0; correctedHR = 0; else; correctedHR = hitRate; end 
                                hitMat.(fld1)(FO,M,a) = correctedHR;
                                out5rows = sig(out5Ix,:);
                                fa = ~isnan(out5rows);
                                anyFA = any(fa,1);
                                gte1FArate = sum(anyFA)/ssSZ;% fixed 6/26/22
                                srchrows = sig(f_ix,:);
                                pts = ~isnan(srchrows);
                                anyPts = any(pts,1);
                                ptsRate = sum(anyPts)/ssSZ;
                                if m == 0; correctedFA = ptsRate; else; correctedFA = gte1FArate; end
                                faMat.(fld1)(FO,M,a) = correctedFA;
                            end %a    
                        end %f1
                        k = k+1; 
						row = cell(1,nV);
						row(1:6) = {s nw f_osc simOpt.pbase*1000-f_osc simOpt.pbase*1000 m};
						mats = {hitMat;faMat};
						col = nnz(~isemptyCell(row))+1;
						for mt = 1:numel(mats)
							for f1 = 1:nF1
								fld1 = fld1arr{f1};
								Mat = mats{mt};
								row(col) = {Mat.(fld1)(FO,M,aIX)};
								col = col+1;
							end %f1
						end %mt						
						for mt = 1:numel(mats)
							Mat = mats{mt};
							for fp = 1:nCO
								fldpair = fldpairs(fp,:);
								fldA = fldpair{1}; fldB = fldpair{2};
								row(col) = {Mat.(fldA)(FO,M,aIX) - Mat.(fldB)(FO,M,aIX)};
								col = col+1;
							end    
						end %mt									
                        T(k,:) = row;
                    end %m
                end %pb
                for f1 = 1:nF1 
                    fld1 = fld1arr{f1};
                    for a = 1:nA
                        hrAll.(fld1)(FN,a,s) = mean(hitMat.(fld1)(:,:,a),'all');
                        faAll.(fld1)(FN,a,s) = mean(faMat.(fld1)(:,:,a),'all');
                    end %a
                end %f1   
            end %ss
        end %FO
end %NW    
toc;
%%
outfile = sprintf('subsample%u_%s.mat',nSS,string(getDateChar('MMddyy')));
outpath = fullfile(datadir,outfile);
fprintf('saving %s (no AUC yet)\n',outpath);
save(outpath)

hrMns = structfun(@(x) squeeze(mean(x,1))',hrAll,'UniformOutput',false);
faMns = structfun(@(x) squeeze(mean(x,1))',faAll,'UniformOutput',false);

% find shared fa subrange
faMaxLB = max(cell2mat(struct2cell(structfun(@(x) max(min(x,[],2)),faMns,'UniformOutput',false))));
faMinUB = min(cell2mat(struct2cell(structfun(@(x) min(max(x,[],2)),faMns,'UniformOutput',false))));

faIX = structfun(@(x) x >= faMaxLB & x <= faMinUB,faMns,'UniformOutput',false);

for f1 = 1:nF1
    fld1 = fld1arr{f1};
    faMatCropFA.(fld1) = cell(nSS,1);
    hrMatCropFA.(fld1) = cell(nSS,1);
    faMatCropHR.(fld1) = cell(nSS,1);
    hrMatCropHR.(fld1) = cell(nSS,1);
    pAUCf.(fld1) = nan(nSS,1);
    pAUCh.(fld1) = nan(nSS,1);
    for s = 1:nSS
        faMn = faMns.(fld1)(s,:);
        hrMn = hrMns.(fld1)(s,:);
        [faMatCropFA.(fld1){s},hrMatCropFA.(fld1){s}] = cropToOneDimRange(faMn,hrMn,faMaxLB,faMinUB);
        [faMatCropHR.(fld1){s},hrMatCropHR.(fld1){s}] = cropToOneDimRange(hrMn,faMn,hrMaxLB,hrMinUB); 
        pAUCf.(fld1)(s) = trapz(faMatCropFA.(fld1){s},hrMatCropFA.(fld1){s});
        pAUCh.(fld1)(s) = trapz(hrMatCropHR.(fld1){s},faMatCropHR.(fld1){s});        
    end %s
end %f
%%
%
for fp = 1:nCO
    fldpair = fldpairs(fp,:);
    fldA = fldpair{1}; fldB = fldpair{2};
    codepair = codepairs(fp,:);
    DiffsF.(codepair) = pAUCf.(fldA)-pAUCf.(fldB);
    mnDiffF.(codepair) = mean(DiffsF.(codepair));
    seDiffF.(codepair) = stderr(DiffsF.(codepair));
    [hF.(codepair),pF.(codepair),ciF.(codepair),statsF.(codepair)]...
        = ttest(pAUCf.(fldA),pAUCf.(fldB));
end %fp




%%
write = true;
outfile = sprintf('subsample%u_%s.mat',nSS,suffix);
outpath = fullfile(datadir,outfile);
if write
	fprintf('saving %s\n',outpath);
	save(outpath); 
end
