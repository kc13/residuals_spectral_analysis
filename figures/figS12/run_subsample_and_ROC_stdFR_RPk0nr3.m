%run_subsample_and_ROC_stdFR_RPk0nr3.m
% assumes run_simulations_stdFR_RPk0nr3.m has been run first

%%
clear all;
close all;
clc;

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
%% 
randdir = pwd;
randfile = 'rnd_info_subsample_stdFR_RPk0nr3.mat';
randpath = fullfile(randdir,randfile);

%% set write prefs down at bottom
outdir = pwd;
suffix = 'stdFR_RPk0nr3_new'; 

%% rnd seed: due to use of parfor, some differences
% in output will remain relative to stored output (see readme)
reuseRndSeed = false;

if reuseRndSeed
	load(randpath)
	rng(rsd,'twister')
end

%% some other relevant vars will load with data files
simOpt.win_len_msec = 1024;
nwArr = [30;60;120];
nNW = numel(nwArr);
foArr = [7 9 12 20 32]';
nFO = numel(foArr);
nFN = nNW*nFO;
fld1arr = {'comp','res'}; 
fldpairs = nchoosek(fld1arr,2);
nF1 = numel(fld1arr);
nIter = 100; % in original simulations: used for constructing filename
codes = ['C','R']';
codepairs = nchoosek(codes,2);
nCO = size(codepairs,1);

pbase_offset = [2.^(0:5)];
nPB = numel(pbase_offset);
m_arr = [0:0.2:1]; 
nM = numel(m_arr);

% for stats
alphaMat = [0.01 0.05]'*(10.^[-6:1:2]);
alphaArr = reshape(alphaMat,[],1);
alphaArr(end) = [];
nA = numel(alphaArr);
cntl_bnds = [250 500];
srch_bnds = [0 100];
fvec = [0:500/512:500]';
f_ix = fvec > srch_bnds(1) & fvec <= srch_bnds(2);
targA = 0.05;
aIX = find(alphaArr == targA);

%% for subsamples
nSS = 1000;
ssSZ = 20; % this has to be downsized also to do 2 iter test
ssIXmat = cell2mat(arrayfun(@(x) randsample(nIter,ssSZ)',[1:nSS]','UniformOutput',false));

%% table preallocate
tvars = [{'SS','NW','FO','PBO','PB','M'},...
    join([repmat({'HR'},nF1,1) cellstr(codes)],'')' join([repmat({'FA'},nF1,1) cellstr(codes)],'')'...
    join([repmat({'HRD'},nCO,1) cellstr(codepairs)],'')' join([repmat({'FAD'},nCO,1) cellstr(codepairs)],'')'];
nV = numel(tvars);
nR = nNW*nFO*nPB*nM*nSS;
vartypes = repmat({'double'},1,nV);
T = table('Size',[nR nV],'VariableNames',tvars,'VariableTypes',vartypes);

k = 0;
%%
dispmod = 100;
tic;
for NW = 1:nNW
    nw = nwArr(NW);
    simOpt.len = nw*simOpt.win_len_msec; 
	for FO = 1:nFO
		FN = (NW-1)*nFO + FO;
		simOpt.f_osc = foArr(FO);
		f_osc = simOpt.f_osc;   
        targIx = knnsearch(fvec,f_osc,'K',3)';
        out5Ix = find(abs(fvec-f_osc)>5);		
		pbase_arr = (simOpt.f_osc + pbase_offset)/1000;
        fname = sprintf('len%uf_osc%uiter%u.mat',simOpt.len,simOpt.f_osc,nIter);
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
			end %f1
			for pb = 1:nPB
				simOpt.pbase = pbase_arr(pb);
				fldstrPB = sprintf('pb%g',simOpt.pbase*1000);
				for M = 1:nM
					m = m_arr(M);
					fldstrM = sprintf('m%g',m*100); 
					freq = sim_data.(fldstrPB).(fldstrM).f;
					for f1 = 1:nF1
						fld1 = fld1arr{f1};
                        S1 = sim_data.(fldstrPB).(fldstrM).(fld1).S1;
                        sig = nan(numel(freq),ssSZ);p_osc = sim_data.(fldstrPB).(fldstrM).p_osc;
						for a = 1:nA
							Alpha = alphaArr(a);
							for si = 1:ssSZ 
                                i = ssIX(si);
                                pow = S1(:,i);
								stats = get_spec_stats(pow,freq,srch_bnds,cntl_bnds,Alpha);sig(:,si) = stats.sigvec1TL;				
							end %s1
							sigrows = sig(targIx,:);
							hits = ~isnan(sigrows);
                            anyHits = any(hits,1);
                            hitRate = sum(anyHits)/ssSZ;      
							if m == 0; correctedHR = 0; else; correctedHR = hitRate; end
							hitMat.(fld1)(pb,M,a) = correctedHR;out5rows = sig(out5Ix,:);
							fa = ~isnan(out5rows);
							anyFA = any(fa,1);						
							gte1FArate = sum(anyFA)/ssSZ;
							srchrows = sig(f_ix,:); 
							pts = ~isnan(srchrows);
							anyPts = any(pts,1);
							ptsRate = sum(anyPts)/ssSZ;
							if m == 0; correctedFA = ptsRate; else; correctedFA = gte1FArate; end
							faMat.(fld1)(pb,M,a) = correctedFA;
						end %a
					end %f1	
					k = k+1;
					row = cell(1,nV);
					row(1:6) = {s nw f_osc simOpt.pbase*1000-f_osc simOpt.pbase*1000 m};
					mats = {hitMat;faMat};
					col = nnz(~isemptyCell(row))+1;
					% columns for individual method results
					for mt = 1:numel(mats)
						for f1 = 1:nF1
                            fld1 = fld1arr{f1};
                            Mat = mats{mt};
                            row(col) = {Mat.(fld1)(pb,M,aIX)};
                            col = col+1;						
						end %f1
					end %mt
					% columns for method differences
                    for mt = 1:numel(mats)
                        Mat = mats{mt};
                        for fp = 1:nCO
                            fldpair = fldpairs(fp,:);
                            fldA = fldpair{1}; fldB = fldpair{2};
                            row(col) = {Mat.(fldA)(pb,M,aIX) - Mat.(fldB)(pb,M,aIX)};
                            col = col+1;
                        end    
                    end %mt					
					T(k,:) = row;				
				end %M				
			end %pb
			for f1 = 1:nF1
				fld1 = fld1arr{f1};
				for a = 1:nA  
					hrAll.(fld1)(FN,a,s) = mean(hitMat.(fld1)(:,:,a),'all');
					faAll.(fld1)(FN,a,s) = mean(faMat.(fld1)(:,:,a),'all'); 
				end %a
			end % f1
		end %s		
	end %FO
end %NW
toc;

% pROC steps below
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
	pAUCf.(fld1) = nan(nSS,1);
	for s = 1:nSS
		faMn = faMns.(fld1)(s,:);
		hrMn = hrMns.(fld1)(s,:);
        [faMatCropFA.(fld1){s},hrMatCropFA.(fld1){s}] = cropToOneDimRange(faMn,hrMn,faMaxLB,faMinUB);		
        pAUCf.(fld1)(s) = trapz(faMatCropFA.(fld1){s},hrMatCropFA.(fld1){s});		
	end %s
end %f1

% t-test comparing method pairs
% written to make general to arbitrary # of methods
% in current practice only one pair
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


write = false;
outfile = sprintf('subsample%u_%s.mat',nSS,suffix);
outpath = fullfile(datadir,outfile);
fprintf('saving %s\n',outpath);
if write; save(outpath); end
